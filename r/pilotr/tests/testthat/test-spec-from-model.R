# Tests for spec_from_model(), which reads a design specification off a fitted lmer. The
# round-trip is the property that matters: a specification simulated, fitted and read back
# has to land on the same design, with estimates near the generating values, because the
# whole point of the function is that the recovered specification can be scaled up and
# trusted to describe the pilot it came from. Kept modest: two small fits, no replicates.

# A small crossed design with a within factor, a between factor and random slopes, and the
# maximal fit of its own modelling data.
sfm_fit <- function() {
  spec <- build_spec(list(name = "pilot", seed = 42, design_kind = "within",
                          include_items = TRUE, n_subject = 30, n_item = 20,
                          factor_name = "cond", lev1 = "a", lev2 = "b",
                          intercept = 6, effect = 0.05,
                          subj_int_sd = 0.12, subj_slope_sd = 0.04, subj_corr = 0.2,
                          item_int_sd = 0.08, item_slope_sd = 0, item_corr = 0,
                          family = "gaussian", resp_name = "", sigma = 0.3))
  spec$factors[[2]] <- list(name = "grp", levels = c("x", "y"),
                            contrasts = list(grp = c(-0.5, 0.5)), between = "subject")
  spec$fixed$coefficients$grp <- 0.1
  d <- simulate_design(spec)
  list(spec = spec,
       fit = lme4::lmer(model_formula(spec), data = model_data(spec, d)))
}

test_that("spec_from_model round-trips units, placement and estimates from a pilot fit", {
  skip_if_not_installed("lme4")
  f <- sfm_fit()
  rec <- suppressMessages(spec_from_model(f$fit, n_subject = 60, n_item = 40))

  # The requested sizes land in the units; without them the fit's own level counts do.
  expect_identical(rec$units$subject$n, 60L)
  expect_identical(rec$units$item$n, 40L)
  asis <- suppressMessages(spec_from_model(f$fit))
  expect_identical(asis$units$subject$n, 30L)
  expect_identical(asis$units$item$n, 20L)

  # Placement: the within factor's contrast column ("effect", build_spec's key) varies
  # inside both units, the between factor's is constant within subjects. Both come back
  # as two-valued numeric columns.
  fac_of <- function(s, contrast) {
    hit <- Filter(function(fc) identical(names(fc$contrasts), contrast), s$factors)
    expect_length(hit, 1L)
    hit[[1]]
  }
  expect_setequal(fac_of(rec, "effect")$vary_within, c("subject", "item"))
  expect_identical(fac_of(rec, "grp")$between, "subject")
  expect_identical(attr(rec, "column_kinds")[["effect"]], "factor")
  expect_identical(attr(rec, "column_kinds")[["grp"]], "factor")

  # The estimates sit near the generating values. The tolerances are wide because one
  # 30-by-20 pilot estimates a variance component with few effective observations, but
  # they still separate the right reading from a wrong field or a wrong scale.
  expect_identical(rec$response$family, "gaussian")
  expect_equal(rec$response$sigma, 0.3, tolerance = 0.1)
  expect_equal(rec$fixed$intercept, 6, tolerance = 0.02)
  expect_equal(rec$fixed$coefficients$effect, 0.05, tolerance = 0.6)
  expect_equal(rec$fixed$coefficients$grp, 0.1, tolerance = 0.6)
  expect_equal(rec$random$subject$intercept_sd, 0.12, tolerance = 0.5)
  expect_equal(rec$random$item$intercept_sd, 0.08, tolerance = 0.5)
  expect_equal(rec$random$subject$slopes$effect, 0.04, tolerance = 0.8)

  # The recovered specification is itself usable: validated on return, simulable now.
  d2 <- simulate_design(rec)
  expect_identical(nrow(d2), 60L * 40L * 2L)   # subjects x items x within-factor levels
})

test_that("spec_from_model re-keys a product column to its interaction", {
  skip_if_not_installed("lme4")
  spec <- load_spec(pilotr_example("reading_time_continuous"))
  spec$units$subject$n <- 24
  spec$units$item$n <- 16
  d <- simulate_design(spec)
  fit <- lme4::lmer(model_formula(spec), data = model_data(spec, d))
  rec <- suppressMessages(spec_from_model(fit, family = "lognormal"))

  # model_data() wrote the interaction "SyntaxPC:age" as the product column "SyntaxPC_age";
  # reading it back as an independent predictor would add a term the design never had.
  keys <- names(rec$fixed$coefficients)
  expect_true("SyntaxPC:age" %in% keys)
  expect_false("SyntaxPC_age" %in% keys)
  expect_false("SyntaxPC_age" %in% names(attr(rec, "column_kinds")))

  # The continuous predictors keep their unit of variation.
  vb <- vapply(rec$predictors, function(p) p$varies_by, character(1))
  names(vb) <- vapply(rec$predictors, function(p) p$name, character(1))
  expect_identical(unname(vb[c("SyntaxPC", "age")]), c("item", "subject"))
  expect_identical(attr(rec, "column_kinds")[["age"]], "predictor")
  expect_identical(rec$response$family, "lognormal")
})

test_that("spec_from_model combines factors between one unit into one factor of cells", {
  skip_if_not_installed("lme4")
  # A balanced 2 x 2 between-subjects pilot, 10 subjects per cell, simulated from the encoding
  # that keeps every cell. Read back as two factors between subjects, it used to give cells of
  # 40, 0, 0 and 40 at 80 subjects, and a refit that dropped B and the interaction.
  pilot <- list(name = "pilot", seed = 11,
                units = list(subject = list(n = 40), item = list(n = 6)),
                factors = list(list(name = "cell", levels = c("a1.b1", "a1.b2", "a2.b1", "a2.b2"),
                                    contrasts = list(a = c(-0.5, -0.5, 0.5, 0.5),
                                                     b = c(-0.5, 0.5, -0.5, 0.5)),
                                    between = "subject")),
                fixed = list(intercept = 10, coefficients = list(a = 0.5, b = 0.3, "a:b" = 0.2)),
                random = list(subject = list(intercept_sd = 1), item = list(intercept_sd = 0.5)),
                response = list(family = "gaussian", name = "y", sigma = 1))
  d <- simulate_design(pilot)
  d$A <- factor(substr(d$cell, 1, 2))
  d$B <- factor(substr(d$cell, 4, 5))
  fit <- lme4::lmer(y ~ A * B + (1 | subject) + (1 | item), data = d)

  said <- character(0)
  keep <- function(m) { said <<- c(said, conditionMessage(m)); invokeRestart("muffleMessage") }
  rec <- withCallingHandlers(spec_from_model(fit, n_subject = 80), message = keep)
  expect_true(any(grepl(paste0(
    "spec_from_model() combined the factors 'A' and 'B', both constant within subject, into one ",
    "factor 'A_B' whose levels are their cells, because pilotr would otherwise assign them to the ",
    "same or overlapping blocks of subjects and confound their effects."), said, fixed = TRUE)))

  # One factor of cells, the first-listed factor varying slowest, carrying both factors' fitted
  # contrast columns, so that every coefficient key the fit produced still names a column.
  expect_length(rec$factors, 1L)
  f <- rec$factors[[1]]
  expect_identical(f$name, "A_B")
  expect_identical(f$between, "subject")
  expect_identical(f$levels, c("a1.b1", "a1.b2", "a2.b1", "a2.b2"))
  expect_identical(f$contrasts, list(Aa2 = c(0, 0, 1, 1), Bb2 = c(0, 1, 0, 1)))
  expect_setequal(names(rec$fixed$coefficients), c("Aa2", "Bb2", "Aa2:Bb2"))

  d2 <- simulate_design(rec)
  expect_identical(as.vector(table(d2$A_B[!duplicated(d2$subject)])), rep(20L, 4))
  said <- character(0)
  refit <- withCallingHandlers(lme4::lmer(model_formula(rec), data = model_data(rec, d2)),
                               message = keep)
  expect_false(any(grepl("rank deficient", said, fixed = TRUE)))
  fe <- lme4::fixef(refit)
  expect_identical(names(fe), c("(Intercept)", "Aa2", "Bb2", "Aa2_Bb2"))
  expect_false(anyNA(fe))
})

test_that("factors between items combine too, keeping each cell's contrasts and a free name", {
  facs <- list(
    list(name = "A", levels = c("a1", "a2"), contrasts = list(Aa2 = c(0, 1)), between = "item"),
    list(name = "cond", levels = c("x", "y"), contrasts = list(cond = c(-0.5, 0.5)),
         vary_within = c("subject", "item")),
    list(name = "B", levels = c("b1", "b2", "b3"),
         contrasts = list(Bb2 = c(0, 1, 0), Bb3 = c(0, 0, 1)), between = "item"))
  expect_message(out <- .combine_between(facs, used = c("A", "B", "cond", "A_B")),
                 paste0("combined the factors 'A' and 'B', both constant within item, into one ",
                        "factor 'A_B_'"), fixed = TRUE)
  # The combined factor takes the first component's place, and the within factor is untouched.
  expect_identical(vapply(out, function(f) f$name, character(1)), c("A_B_", "cond"))
  expect_identical(out[[2]], facs[[2]])
  f <- out[[1]]
  expect_identical(f$between, "item")
  expect_identical(f$levels, c("a1.b1", "a1.b2", "a1.b3", "a2.b1", "a2.b2", "a2.b3"))
  expect_identical(f$contrasts, list(Aa2 = c(0, 0, 0, 1, 1, 1), Bb2 = c(0, 1, 0, 0, 1, 0),
                                     Bb3 = c(0, 0, 1, 0, 0, 1)))
})

test_that("spec_from_model refuses what it cannot read, saying what to do instead", {
  skip_if_not_installed("lme4")
  # No random effects: the part of a specification hardest to guess is missing.
  lmfit <- stats::lm(y ~ x, data = data.frame(x = 1:20, y = rnorm(20)))
  expect_error(spec_from_model(lmfit), "A model with no random effects", fixed = TRUE)
  # A generalised model keeps its random effects on the link scale.
  gdat <- data.frame(y = rep(c(0L, 1L), 20), g = factor(rep(1:8, each = 5)))
  gfit <- suppressWarnings(suppressMessages(
    lme4::glmer(y ~ 1 + (1 | g), data = gdat, family = stats::binomial())))
  expect_error(spec_from_model(gfit), "fitted with glmer()", fixed = TRUE)
  # Not a model at all.
  expect_error(spec_from_model(42), "not an object of class 'numeric'", fixed = TRUE)
})
