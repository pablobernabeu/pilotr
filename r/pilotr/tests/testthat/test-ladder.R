# A preregistered simplification rule. power_mixed() and precision_design() fitted one model to
# every replicate, by default the maximal model, which is singular in most replicates of a small
# crossed design. An analysis plan that drops the correlations, then the slopes, whenever the fit is
# singular could not be simulated, although the rule changes power and the Type I error rate. A
# ladder of formulas, with a fallback rule, simulates it: each replicate keeps the first rung whose
# fit the rule accepts, or else the last, and the result reports the share of replicates analysed
# at each rung.

# The crossed priming design of the vignettes, shrunk to 8 subjects and 6 items, where the maximal
# model is singular in most replicates.
priming <- function(n_subject = 8, n_item = 6) build_spec(list(
  name = "priming", seed = 1, design_kind = "within", include_items = TRUE,
  n_subject = n_subject, n_item = n_item,
  factor_name = "condition", lev1 = "related", lev2 = "unrelated",
  intercept = 6, effect = 0.05,
  subj_int_sd = 0.12, subj_slope_sd = 0.04, subj_corr = 0.2,
  item_int_sd = 0.08, item_slope_sd = 0.02, item_corr = -0.1,
  family = "shifted_lognormal", resp_name = "RT", sigma = 0.3, shift = 200))

LADDER <- c("maximal", "zcp", "intercepts")

test_that("model_formula() builds the maximal, zero-correlation and intercepts-only models", {
  spec <- priming()
  expect_identical(model_formula(spec), model_formula(spec, "maximal"))
  expect_identical(deparse(model_formula(spec)),
                   ".y ~ effect + (1 + effect | subject) + (1 + effect | item)")
  expect_identical(deparse(model_formula(spec, "zcp")),
                   ".y ~ effect + (1 + effect || subject) + (1 + effect || item)")
  expect_identical(deparse(model_formula(spec, "intercepts")),
                   ".y ~ effect + (1 | subject) + (1 | item)")
  # A grouping factor without slopes has an intercept alone in every structure, and a
  # specification that already leaves its slopes uncorrelated has a maximal model with `||`.
  spec$random$item$slopes <- NULL
  spec$random$item$correlations <- NULL
  spec$random$subject$correlations <- NULL
  expect_identical(deparse(model_formula(spec)),
                   ".y ~ effect + (1 + effect || subject) + (1 | item)")
  expect_identical(model_formula(spec, "zcp"), model_formula(spec))
  expect_error(model_formula(spec, "minimal"), "should be one of")
})

test_that("each replicate keeps the first rung whose fit the rule accepts, else the last", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  spec <- priming()
  ladder <- pilotr:::.as_ladder(LADDER, "singular", spec)
  seeds <- replicate_seeds(spec$seed, 6)
  prep <- function(d) model_data(spec, d)
  kept <- character(0)
  for (i in seq_along(seeds)) {
    r <- pilotr:::.design_rep(i, spec, seeds, prep, ladder, "effect")
    s <- spec
    s$seed <- seeds[i]
    d <- prep(simulate_design(s))
    fits <- lapply(LADDER, function(k) pilotr:::.fit_model(model_formula(spec, k), d,
                                                           family = "shifted_lognormal",
                                                           test = TRUE, auto = TRUE))
    ok <- vapply(fits, function(f) !is.null(f$fit) && !f$singular, logical(1))
    k <- if (any(ok)) which(ok)[1] else length(LADDER)
    expect_identical(r$rung, LADDER[k], info = i)
    expect_identical(r$est[["effect"]], lme4::fixef(fits[[k]]$fit)[["effect"]], info = i)
    expect_identical(r$singular, fits[[k]]$singular, info = i)
    kept <- c(kept, r$rung)
  }
  # The maximal model is singular often enough here that the rule moves down the ladder.
  expect_true(any(kept != "maximal"))
})

test_that("fallback = 'warning' also moves on from a fit with a warning", {
  passes <- pilotr:::.rung_passes
  clean <- list(fit = TRUE, singular = FALSE, messages = character(0), converged = TRUE)
  warned <- list(fit = TRUE, singular = FALSE, messages = "Model failed to converge",
                 converged = FALSE)
  singular <- list(fit = TRUE, singular = TRUE, messages = character(0), converged = FALSE)
  failed <- list(fit = NULL, singular = FALSE, messages = "an error", converged = FALSE)
  expect_true(passes(clean, "singular"))
  expect_true(passes(clean, "warning"))
  expect_true(passes(warned, "singular"))
  expect_false(passes(warned, "warning"))
  expect_false(passes(singular, "singular"))
  expect_false(passes(singular, "warning"))
  expect_false(passes(failed, "singular"))
  expect_false(passes(failed, "warning"))
})

test_that("a singular-prone design's ladder reaches later rungs, with shares that sum to 1", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  spec <- priming()
  maximal <- power_mixed(spec, n_sims = 8)
  ladder <- power_mixed(spec, formula = LADDER, n_sims = 8)
  expect_identical(names(ladder$rung_shares), LADDER)
  expect_equal(sum(ladder$rung_shares), 1)
  expect_gt(sum(ladder$rung_shares[c("zcp", "intercepts")]), 0)
  expect_identical(ladder$fallback, "singular")
  # The same replicates, so the maximal rung keeps exactly the fits that were not singular.
  expect_equal(ladder$rung_shares[["maximal"]],
               (maximal$n_returned - maximal$n_singular) / 8)
  # The counts describe the fit each replicate kept, which is singular only at the last rung.
  expect_lte(ladder$n_singular, maximal$n_singular)
  expect_identical(ladder$fitter, maximal$fitter)

  pr <- precision_design(spec, formula = LADDER, rope = 0.02, n_sims = 8)
  shares <- unlist(pr[1, paste0("rung_", LADDER)])
  expect_equal(sum(shares), 1)
  expect_gt(sum(shares[c("rung_zcp", "rung_intercepts")]), 0)
})

test_that("a single formula returns what it returned before", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  spec <- priming()
  single <- power_mixed(spec, n_sims = 3)
  expect_null(single$rung_shares)
  expect_null(single$fallback)
  # The rule applies only to a ladder.
  expect_identical(power_mixed(spec, n_sims = 3, fallback = "warning"), single)
  # A ladder of one rung is that model's analysis, with the share it adds.
  one <- power_mixed(spec, formula = "maximal", n_sims = 3)
  expect_identical(one$rung_shares, c(maximal = 1))
  expect_identical(unclass(one)[names(single)], unclass(single))
  pr <- precision_design(spec, rope = 0.02, n_sims = 3)
  expect_false(any(startsWith(names(pr), "rung_")))
})

test_that("rungs named as structures take the family's model, and formulas are fitted as written", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  spec <- build_spec(list(
    name = "accuracy", seed = 1, design_kind = "within", include_items = TRUE,
    n_subject = 10, n_item = 8, factor_name = "cond", lev1 = "a", lev2 = "b",
    intercept = 1, effect = 0.5, subj_int_sd = 0.5, subj_slope_sd = 0.2, subj_corr = 0.2,
    item_int_sd = 0.3, item_slope_sd = 0.1, item_corr = 0, family = "bernoulli",
    resp_name = ""))
  derived <- suppressWarnings(power_mixed(spec, formula = c("maximal", "intercepts"), n_sims = 2))
  expect_identical(derived$fitter, "lme4::glmer (binomial), Wald z")
  given <- suppressWarnings(power_mixed(
    spec, formula = list(maximal = model_formula(spec),
                         intercepts = model_formula(spec, "intercepts")), n_sims = 2))
  expect_identical(given$fitter, "lmerTest::lmer, Satterthwaite t")
  expect_identical(names(given$rung_shares), c("maximal", "intercepts"))
  # Rungs whose fitters differ are each named.
  mixed <- power_mixed(spec, formula = list(.y ~ effect + (1 | subject), .y ~ effect), n_sims = 2)
  expect_identical(mixed$fitter, "lmerTest::lmer, Satterthwaite t; stats::lm, t")
  expect_identical(names(mixed$rung_shares), c("1", "2"))
})

test_that("a one-sided test applies to the fit each replicate kept", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  spec <- priming()
  ladder <- pilotr:::.as_ladder(LADDER, "singular", spec)
  seeds <- replicate_seeds(spec$seed, 4)
  prep <- function(d) model_data(spec, d)
  at <- function(i, alternative)
    pilotr:::.design_rep(i, spec, seeds, prep, ladder, "effect", alternative = alternative)
  for (i in seq_along(seeds)) {
    two <- at(i, "two.sided")
    toward <- if (two$est[["effect"]] > 0) "greater" else "less"
    one <- at(i, toward)
    # The rule chooses the rung before any test, so the direction leaves it and the estimate alone.
    expect_identical(one$rung, two$rung, info = i)
    expect_identical(one$est, two$est, info = i)
    expect_equal(one$p[["effect"]], two$p[["effect"]] / 2, info = i)
  }
  # Every rung's fitter is named with the direction of the test.
  mixed <- power_mixed(spec, formula = list(.y ~ effect + (1 | subject), .y ~ effect), n_sims = 2,
                       alternative = "greater")
  expect_identical(mixed$fitter, paste("lmerTest::lmer, Satterthwaite t, one-sided (greater);",
                                       "stats::lm, t, one-sided (greater)"))
})

test_that("a sweep carries the share of each rung on every row", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  curve <- sweep_spec(priming(), "units$subject$n", c(8, 10), power_mixed,
                      formula = c("maximal", "intercepts"), n_sims = 2)
  expect_true(all(c("rung_maximal", "rung_intercepts") %in% names(curve)))
  expect_equal(curve$rung_maximal + curve$rung_intercepts, c(1, 1))
  prc <- precision_curve(priming(), subject_ns = c(8, 10), formula = c("maximal", "intercepts"),
                         fallback = "warning", rope = 0.02, n_sims = 2)
  expect_equal(prc$rung_maximal + prc$rung_intercepts, c(1, 1))
  expect_error(precision_curve(priming(), subject_ns = 8, formula = c("maximal", "none")),
               "'none' is not a structure")
  expect_error(precision_curve(priming(), subject_ns = 8, fallback = "none"), "should be one of")
})

test_that("a malformed ladder is refused before any fitting", {
  # Both functions ask for their fitting packages before they read their arguments.
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  spec <- priming()
  expect_error(power_mixed(spec, formula = c("maximal", "minimal")),
               "'minimal' is not a structure")
  expect_error(power_mixed(spec, formula = c("maximal", "maximal")), "more than once")
  expect_error(power_mixed(spec, formula = list(a = .y ~ effect, a = .y ~ 1)), "more than once")
  expect_error(power_mixed(spec, formula = list("maximal", "zcp")), "a formula")
  expect_error(power_mixed(spec, formula = list()), "at least one rung")
  expect_error(precision_design(spec, formula = c("zcp", "none")), "'none' is not a structure")
  expect_error(power_mixed(spec, fallback = "converged"), "should be one of")
  expect_error(precision_design(spec, fallback = "converged"), "should be one of")
})
