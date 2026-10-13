# Tests for spec_from_model(), which reads a design specification off a fitted lmer. The
# round-trip is the property that matters: a specification simulated, fitted and read back
# has to land on the same design, with estimates near the generating values, because the
# whole point of the function is that the recovered specification can be scaled up and
# trusted to describe the pilot it came from. Kept modest: two small fits, no replicates.

# A small crossed design with a within factor, a between factor and random slopes, and the
# maximal fit of its own modelling data. `grp` is the true effect of the between factor, and
# `lmer` the fitter, lme4's or lmerTest's.
sfm_fit <- function(grp = 0.1, lmer = lme4::lmer) {
  spec <- build_spec(list(name = "pilot", seed = 42, design_kind = "within",
                          include_items = TRUE, n_subject = 30, n_item = 20,
                          factor_name = "cond", lev1 = "a", lev2 = "b",
                          intercept = 6, effect = 0.05,
                          subj_int_sd = 0.12, subj_slope_sd = 0.04, subj_corr = 0.2,
                          item_int_sd = 0.08, item_slope_sd = 0, item_corr = 0,
                          family = "gaussian", resp_name = "", sigma = 0.3))
  spec$factors[[2]] <- list(name = "grp", levels = c("x", "y"),
                            contrasts = list(grp = c(-0.5, 0.5)), between = "subject")
  spec$fixed$coefficients$grp <- grp
  d <- simulate_design(spec)
  list(spec = spec,
       fit = lmer(model_formula(spec), data = model_data(spec, d)))
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

# The messages a call raises, muffled and without message()'s closing newline, together with the
# call's value.
sfm_said <- function(expr) {
  said <- character(0)
  value <- withCallingHandlers(expr, message = function(m) {
    said <<- c(said, sub("\n$", "", conditionMessage(m)))
    invokeRestart("muffleMessage")
  })
  list(value = value, said = said)
}

# The fit's estimates and standard errors, keyed by term.
sfm_estimates <- function(fit) {
  list(est = lme4::fixef(fit), se = sqrt(diag(as.matrix(stats::vcov(fit)))))
}

test_that("spec_from_model says which fixed effects it took from the pilot, and keeps them", {
  skip_if_not_installed("lme4")
  f <- sfm_fit()
  e <- sfm_estimates(f$fit)
  out <- sfm_said(spec_from_model(f$fit))
  rec <- out$value

  # The default is unchanged: the coefficients are the pilot's point estimates, as they were.
  expect_identical(rec$fixed$coefficients, list(effect = e$est[["effect"]], grp = e$est[["grp"]]))
  expect_identical(attr(rec, "effects"),
                   list(method = "estimate", level = NA_real_, sesoi = character(0)))
  expect_identical(suppressMessages(spec_from_model(f$fit, effects = "estimate")), rec)

  # One message lists them with their values and recommends a SESOI or a safeguard value.
  said <- grep("as the pilot estimated them", out$said, value = TRUE, fixed = TRUE)
  expect_length(said, 1L)
  expect_match(said, sprintf("took 'effect' = %s and 'grp' = %s from the fit",
                             format(e$est[["effect"]], digits = 3),
                             format(e$est[["grp"]], digits = 3)), fixed = TRUE)
  expect_match(said, "smallest effect size of interest with `sesoi`", fixed = TRUE)
  expect_match(said, "effects = \"safeguard\"", fixed = TRUE)
})

test_that("effects = 'safeguard' takes the bound of each 60% interval nearer zero", {
  skip_if_not_installed("lme4")
  f <- sfm_fit()
  e <- sfm_estimates(f$fit)
  out <- sfm_said(spec_from_model(f$fit, effects = "safeguard"))
  rec <- out$value

  # Both estimates are positive, so the bound nearer zero is the lower one, est - qnorm(0.8) se.
  for (k in c("effect", "grp")) {
    v <- rec$fixed$coefficients[[k]]
    expect_gt(v, 0)
    expect_lt(v, e$est[[k]])
    expect_equal(v, e$est[[k]] - stats::qnorm(0.8) * e$se[[k]], tolerance = 1e-12)
  }
  expect_identical(attr(rec, "effects"),
                   list(method = "safeguard", level = 0.6, sesoi = character(0)))
  said <- grep("safeguard values", out$said, value = TRUE, fixed = TRUE)
  expect_length(said, 1L)
  expect_match(said, sprintf("took 'effect' = %s and 'grp' = %s from the fit as safeguard values",
                             format(rec$fixed$coefficients$effect, digits = 3),
                             format(rec$fixed$coefficients$grp, digits = 3)), fixed = TRUE)
  expect_match(said, "60% interval, computed with the normal distribution", fixed = TRUE)
  expect_false(grepl("tends to be underpowered", said, fixed = TRUE))

  # Nothing else changes: the intercept, the variance components and the design are the fit's.
  asis <- suppressMessages(spec_from_model(f$fit))
  rec$fixed$coefficients <- asis$fixed$coefficients
  attr(rec, "effects") <- attr(asis, "effects") <- NULL
  expect_identical(rec, asis)
  # The level can be changed. At 70%, both bounds move further from the estimates and stay above
  # zero, grp's by about a third of its standard error.
  wide <- suppressMessages(spec_from_model(f$fit, effects = "safeguard", safeguard_level = 0.7))
  for (k in c("effect", "grp"))
    expect_equal(wide$fixed$coefficients[[k]],
                 e$est[[k]] - stats::qnorm(0.85) * e$se[[k]], tolerance = 1e-12)
  expect_identical(attr(wide, "effects")$level, 0.7)
})

test_that("the safeguard value of a negative effect is its upper bound", {
  skip_if_not_installed("lme4")
  f <- sfm_fit(grp = -0.2)
  e <- sfm_estimates(f$fit)
  v <- suppressMessages(spec_from_model(f$fit, effects = "safeguard"))$fixed$coefficients$grp
  expect_lt(e$est[["grp"]], v)
  expect_lt(v, 0)
  expect_equal(v, e$est[["grp"]] + stats::qnorm(0.8) * e$se[["grp"]], tolerance = 1e-12)
})

test_that("a lmerTest fit takes its safeguard values from Satterthwaite's t", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  fit <- sfm_fit(lmer = lmerTest::lmer)$fit
  expect_s4_class(fit, "lmerModLmerTest")
  co <- summary(fit)$coefficients
  out <- sfm_said(spec_from_model(fit, effects = "safeguard"))
  for (k in c("effect", "grp"))
    expect_equal(out$value$fixed$coefficients[[k]],
                 co[k, "Estimate"] - stats::qt(0.8, co[k, "df"]) * co[k, "Std. Error"],
                 tolerance = 1e-12)
  expect_match(grep("safeguard values", out$said, value = TRUE, fixed = TRUE),
               "computed with the t distribution with Satterthwaite's degrees of freedom",
               fixed = TRUE)
})

test_that("a lmerTest fit without Satterthwaite's degrees of freedom takes the normal quantile", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  fit <- sfm_fit(lmer = lmerTest::lmer)$fit
  # lmerTest computes the degrees of freedom from the covariance of the variance parameters that
  # it keeps in the fit. Blanking that matrix leaves summary() reporting them as NA.
  skip_if_not(methods::.hasSlot(fit, "vcov_varpar"), "this lmerTest keeps no vcov_varpar slot")
  fit@vcov_varpar[] <- NA_real_
  e <- sfm_estimates(fit)
  out <- sfm_said(spec_from_model(fit, effects = "safeguard"))
  for (k in c("effect", "grp"))
    expect_equal(out$value$fixed$coefficients[[k]],
                 e$est[[k]] - stats::qnorm(0.8) * e$se[[k]], tolerance = 1e-12)
  expect_match(grep("safeguard values", out$said, value = TRUE, fixed = TRUE),
               "computed with the normal distribution", fixed = TRUE)
})

test_that("an effect whose 60% interval includes zero stops the safeguard, pointing to sesoi", {
  skip_if_not_installed("lme4")
  # With a true effect of 0.05, this pilot estimates grp at about 0.003 with a standard error of
  # about 0.039, so its 60% interval holds zero and the pilot leaves even the sign open. A value of
  # zero would plan a study of no effect at all, so the call stops.
  f <- sfm_fit(grp = 0.05)
  err <- tryCatch(suppressMessages(spec_from_model(f$fit, effects = "safeguard")),
                  error = identity)
  expect_s3_class(err, "error")
  msg <- conditionMessage(err)
  expect_match(msg, "cannot give 'grp' a safeguard value: its 60% interval, from ", fixed = TRUE)
  expect_match(msg, "includes zero", fixed = TRUE)
  expect_match(msg, "sesoi = c(grp = ...)", fixed = TRUE)
  expect_false(grepl("'effect'", msg, fixed = TRUE))

  # A smallest effect size of interest for grp lets the safeguard go ahead for the other effect.
  out <- sfm_said(spec_from_model(f$fit, effects = "safeguard", sesoi = c(grp = 0.04)))
  rec <- out$value
  expect_identical(rec$fixed$coefficients$grp, 0.04)
  e <- sfm_estimates(f$fit)
  expect_equal(rec$fixed$coefficients$effect,
               e$est[["effect"]] - stats::qnorm(0.8) * e$se[["effect"]], tolerance = 1e-12)
  expect_identical(attr(rec, "effects"), list(method = "safeguard", level = 0.6, sesoi = "grp"))
  expect_true(any(grepl("set 'grp' = 0.04 from `sesoi` and took 'effect' = ", out$said,
                        fixed = TRUE)))
})

test_that("sesoi sets the effects it names and leaves the others as estimated", {
  skip_if_not_installed("lme4")
  f <- sfm_fit()
  e <- sfm_estimates(f$fit)
  out <- sfm_said(spec_from_model(f$fit, sesoi = c(grp = 0.03)))
  rec <- out$value
  expect_identical(rec$fixed$coefficients, list(effect = e$est[["effect"]], grp = 0.03))
  expect_identical(attr(rec, "effects"), list(method = "estimate", level = NA_real_, sesoi = "grp"))
  said <- grep("from `sesoi`", out$said, value = TRUE, fixed = TRUE)
  expect_length(said, 1L)
  expect_match(said, sprintf(paste0("set 'grp' = 0.03 from `sesoi` and took 'effect' = %s from ",
                                    "the fit, as the pilot estimated it."),
                             format(e$est[["effect"]], digits = 3)), fixed = TRUE)
  # With every effect named, nothing is left to take from the fit, and nothing to recommend.
  out <- sfm_said(spec_from_model(f$fit, sesoi = c(grp = 0.03, effect = 0.02)))
  expect_identical(out$value$fixed$coefficients, list(effect = 0.02, grp = 0.03))
  said <- grep("from `sesoi`", out$said, value = TRUE, fixed = TRUE)
  expect_identical(said, "spec_from_model() set 'effect' = 0.02 and 'grp' = 0.03 from `sesoi`.")
})

test_that("spec_from_model refuses an unknown sesoi name and malformed effect arguments", {
  skip_if_not_installed("lme4")
  f <- sfm_fit()
  # The names are checked once the fit's terms have their keys, after the messages about how the
  # columns were read.
  expect_error(suppressMessages(spec_from_model(f$fit, sesoi = c(cond = 0.03))),
               paste0("`sesoi` names 'cond', which is not a coefficient of the returned ",
                      "specification. Its coefficients are 'effect' and 'grp'."), fixed = TRUE)
  # The intercept is a location, which a specification keeps out of its coefficients.
  expect_error(suppressMessages(spec_from_model(f$fit, sesoi = c(intercept = 6, cond = 0.03))),
               "`sesoi` names 'intercept' and 'cond', which are not coefficients", fixed = TRUE)
  for (bad in list(0.03, c(grp = NA), c(grp = Inf), list(grp = 0.03), c(grp = "0.03"),
                   c(0.1, grp = 0.03), c(grp = 0.03, grp = 0.05), numeric(0)))
    expect_error(spec_from_model(f$fit, sesoi = bad), "`sesoi` must be NULL or", fixed = TRUE)
  for (bad in list(0, 1, -0.5, 1.5, NA_real_, "0.6", c(0.6, 0.8)))
    expect_error(spec_from_model(f$fit, effects = "safeguard", safeguard_level = bad),
                 "`safeguard_level` must be a single number greater than 0 and less than 1",
                 fixed = TRUE)
  # No option sets the effects to zero, since design_conditions() builds that null design for a
  # sweep.
  expect_error(spec_from_model(f$fit, effects = "zero"), "should be one of")
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

  # `sesoi` is keyed as the returned specification is, so the interaction is "SyntaxPC:age", and
  # the product column's name is refused.
  set <- suppressMessages(spec_from_model(fit, family = "lognormal",
                                          sesoi = c("SyntaxPC:age" = 0.01)))
  expect_identical(set$fixed$coefficients[["SyntaxPC:age"]], 0.01)
  expect_error(suppressMessages(
    spec_from_model(fit, family = "lognormal", sesoi = c(SyntaxPC_age = 0.01))),
    "`sesoi` names 'SyntaxPC_age', which is not a coefficient", fixed = TRUE)
})

test_that("the sesoi example in a refusal quotes a key that is not a syntactic name", {
  expect_identical(.sesoi_example(c("grp", "SyntaxPC:age")),
                   "sesoi = c(grp = ..., `SyntaxPC:age` = ...)")
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
  # The message gives the reason that holds from specification version 0.4, which allocates two
  # between factors jointly: the factor of cells keeps the specification readable by releases
  # that read only 0.3. It used to say that the factors would otherwise confound their effects.
  expect_true(any(grepl(paste0(
    "spec_from_model() combined the factors 'A' and 'B', both constant within subject, into one ",
    "factor 'A_B' whose levels are their cells. It allocates the subjects as the separate factors ",
    "would under specification version 0.4, and releases that read only version 0.3 read it too."),
    said, fixed = TRUE)))

  # One factor of cells, the first-listed factor varying slowest, carrying both factors' fitted
  # contrast columns, so that every coefficient key the fit produced still names a column.
  expect_length(rec$factors, 1L)
  f <- rec$factors[[1]]
  expect_identical(f$name, "A_B")
  expect_identical(f$between, "subject")
  expect_identical(f$levels, c("a1.b1", "a1.b2", "a2.b1", "a2.b2"))
  expect_identical(f$contrasts, list(Aa2 = c(0, 0, 1, 1), Bb2 = c(0, 1, 0, 1)))
  expect_setequal(names(rec$fixed$coefficients), c("Aa2", "Bb2", "Aa2:Bb2"))
  # Every release that reads version 0.3 reads a factor of cells, so the specification declares
  # 0.3, where two factors between subjects would have needed 0.4.
  expect_identical(names(rec)[1], "spec_version")
  expect_identical(rec$spec_version, "0.3")

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
                        "factor 'A_B_' whose levels are their cells. It allocates the items as ",
                        "the separate factors would"), fixed = TRUE)
  # The combined factor takes the first component's place, and the within factor is untouched.
  expect_identical(vapply(out, function(f) f$name, character(1)), c("A_B_", "cond"))
  expect_identical(out[[2]], facs[[2]])
  f <- out[[1]]
  expect_identical(f$between, "item")
  expect_identical(f$levels, c("a1.b1", "a1.b2", "a1.b3", "a2.b1", "a2.b2", "a2.b3"))
  expect_identical(f$contrasts, list(Aa2 = c(0, 0, 0, 1, 1, 1), Bb2 = c(0, 1, 0, 0, 1, 0),
                                     Bb3 = c(0, 0, 1, 0, 0, 1)))
})

# The crossed priming example, simulated in full and then cut down to the rows a pilot ran.
# `keep` takes the data, whose `level` column holds each row's level index (0 or 1), and returns
# the rows to keep. The pilot carries the outcome as the log of the shifted reaction time, which
# is Gaussian under the example's shifted lognormal family. It carries the condition as a
# -0.5/0.5 column, as a user's own pilot data would.
sfm_priming_pilot <- function(keep) {
  spec <- load_spec(pilotr_example("crossed_mixed_rt"))
  d <- simulate_design(spec)
  d$level <- match(d$condition, spec$factors[[1]]$levels) - 1L
  d <- d[keep(d), ]
  d$y <- log(d$RT - spec$response$shift)
  d$cond <- d$level - 0.5
  d
}

test_that("spec_from_model warns that a counterbalanced pilot comes back fully crossed", {
  skip_if_not_installed("lme4")
  # The pilot rotates two lists in a Latin square. Each subject sees each item once, and the
  # level alternates over items within a subject and over subjects within an item, so every
  # item appears in both conditions across subjects. This is how priming studies are run.
  pilot <- sfm_priming_pilot(function(d) ((d$subject %% 2L) + d$item + d$level) %% 2L == 0L)
  expect_identical(nrow(pilot), 30L * 24L)
  fit <- suppressMessages(suppressWarnings(
    lme4::lmer(y ~ cond + (1 + cond | subject) + (1 + cond | item), data = pilot)))
  expect_warning(
    rec <- suppressMessages(spec_from_model(fit, n_subject = 30, n_item = 24)),
    paste0("spec_from_model(): in the pilot each subject saw each item under one level of ",
           "'cond' (a counterbalanced, Latin-square design), but the returned specification ",
           "crosses every subject with every item under every level, doubling the observations ",
           "per subject and overstating power. To plan the pilot's design, set counterbalance to ",
           "true on 'cond' and declare spec_version \"0.4\" (see Counterbalanced lists in ",
           "SPEC.md). For a release that reads no version above 0.3, encode it as list (between ",
           "subject) x item set (between item), as in SPEC.md's two-list example."),
    fixed = TRUE)
  # As the warning says, the specification crosses every pair with both levels, so it simulates
  # twice the pilot's rows. Changed as the warning advises, it simulates as many as the pilot had.
  expect_setequal(rec$factors[[1]]$vary_within, c("subject", "item"))
  expect_identical(nrow(simulate_design(rec)), 2L * nrow(pilot))
  rec$spec_version <- "0.4"
  rec$factors[[1]]$counterbalance <- TRUE
  expect_identical(nrow(simulate_design(rec)), nrow(pilot))
})

test_that("a crossed pilot with missing rows, or a factor between items, draws no such warning", {
  skip_if_not_installed("lme4")
  # Every tenth row dropped from the full crossing leaves most subject-item pairs with both
  # levels, so the pilot was crossed, however incomplete.
  crossed <- sfm_priming_pilot(function(d) seq_len(nrow(d)) %% 10L != 0L)
  fit <- suppressMessages(suppressWarnings(
    lme4::lmer(y ~ cond + (1 + cond | subject) + (1 + cond | item), data = crossed)))
  expect_no_warning(suppressMessages(spec_from_model(fit)))
  # Here each item has one level for every subject, so the items carry the condition and the
  # factor is placed between items.
  by_item <- sfm_priming_pilot(function(d) (d$item + d$level) %% 2L == 0L)
  fit <- suppressMessages(suppressWarnings(
    lme4::lmer(y ~ cond + (1 + cond | subject) + (1 | item), data = by_item)))
  expect_no_warning(rec <- suppressMessages(spec_from_model(fit)))
  expect_identical(rec$factors[[1]]$between, "item")
})

test_that("the counterbalancing warning counts the levels of a factor with more than two", {
  skip_if_not_installed("lme4")
  # Three primes rotate over three lists and are read from a character column. The
  # specification crosses each subject-item pair with all three, three times the pilot's
  # observations.
  three <- list(name = "three", seed = 5,
                units = list(subject = list(n = 18), item = list(n = 12)),
                factors = list(list(name = "prime", levels = c("related", "neutral", "unrelated"),
                                    contrasts = list(p2 = c(0, 1, 0), p3 = c(0, 0, 1)),
                                    vary_within = c("subject", "item"))),
                fixed = list(intercept = 6, coefficients = list(p2 = 0.03, p3 = 0.06)),
                random = list(subject = list(intercept_sd = 0.1), item = list(intercept_sd = 0.05)),
                response = list(family = "gaussian", name = "y", sigma = 0.3))
  d <- simulate_design(three)
  level <- match(d$prime, three$factors[[1]]$levels) - 1L
  pilot <- d[(d$subject + d$item + level) %% 3L == 0L, ]
  fit <- suppressMessages(suppressWarnings(
    lme4::lmer(y ~ prime + (1 | subject) + (1 | item), data = pilot)))
  w <- expect_warning(rec <- suppressMessages(spec_from_model(fit)), "one level of 'prime'",
                      fixed = TRUE)
  expect_match(conditionMessage(w),
               "every level, multiplying the observations per subject by 3 and overstating power",
               fixed = TRUE)
  # The two-list encoding writes two levels only, so the warning names the field alone.
  expect_true(endsWith(conditionMessage(w), paste0(
    "To plan the pilot's design, set counterbalance to true on 'prime' and declare spec_version ",
    "\"0.4\" (see Counterbalanced lists in SPEC.md).")))
  rec$spec_version <- "0.4"
  rec$factors[[1]]$counterbalance <- TRUE
  expect_identical(nrow(simulate_design(rec)), nrow(pilot))
})

test_that("factors counterbalanced together draw one warning with their joint multiple", {
  skip_if_not_installed("lme4")
  # A 2 x 2 design rotated over four lists: each subject sees each item in one of the four cells.
  # The specification crosses every pair with all four cells, four times the pilot's rows, which
  # one warning per factor, each saying "doubling", would understate.
  cells <- list(name = "cells", seed = 3,
                units = list(subject = list(n = 16), item = list(n = 16)),
                factors = list(
                  list(name = "A", levels = c("a1", "a2"), contrasts = list(a = c(-0.5, 0.5)),
                       vary_within = c("subject", "item")),
                  list(name = "B", levels = c("b1", "b2"), contrasts = list(b = c(-0.5, 0.5)),
                       vary_within = c("subject", "item"))),
                fixed = list(intercept = 6, coefficients = list(a = 0.05, b = 0.05)),
                random = list(subject = list(intercept_sd = 0.1), item = list(intercept_sd = 0.05)),
                response = list(family = "gaussian", name = "y", sigma = 0.3))
  d <- simulate_design(cells)
  cell <- 2L * (match(d$A, c("a1", "a2")) - 1L) + match(d$B, c("b1", "b2")) - 1L
  pilot <- d[(d$subject + d$item + cell) %% 4L == 0L, ]
  fit <- suppressMessages(suppressWarnings(
    lme4::lmer(y ~ A * B + (1 | subject) + (1 | item), data = pilot)))
  w <- capture_warnings(rec <- suppressMessages(spec_from_model(fit)))
  expect_length(w, 1L)
  expect_match(w, paste0("each item under one combination of the levels of 'A' and 'B' (a ",
                         "counterbalanced, Latin-square design), but the returned specification ",
                         "crosses every subject with every item under every combination, ",
                         "multiplying the observations per subject by 4 and overstating power. ",
                         "To plan the pilot's design, set counterbalance to true on 'A' and 'B' ",
                         "and declare spec_version \"0.4\" (see Counterbalanced lists in ",
                         "SPEC.md)."),
               fixed = TRUE)
  expect_false(grepl("two-list", w, fixed = TRUE))
  expect_identical(nrow(simulate_design(rec)), 4L * nrow(pilot))
  # Counterbalanced together, the two factors rotate over the four lists, one cell per pair.
  rec$spec_version <- "0.4"
  expect_identical(vapply(rec$factors, `[[`, "", "name"), c("A", "B"))
  for (k in 1:2) rec$factors[[k]]$counterbalance <- TRUE
  expect_identical(nrow(simulate_design(rec)), nrow(pilot))
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
