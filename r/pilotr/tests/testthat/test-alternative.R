# A directional prediction is tested one-sided, and `alternative` plans that test. Each power
# function took only a two-sided test. The one-sided p-value comes from the statistic behind the
# two-sided one, so it is half the two-sided p-value when the estimate lies in the tested
# direction. On the same replicates, a one-sided test at alpha therefore rejects exactly where the
# two-sided test at 2 alpha rejects with the estimate in that direction. One-sided power at alpha
# is then two-sided power at 2 alpha times one minus Type S at 2 alpha, which the tests below hold
# in counts, exactly.

# A two-group design with a small effect, so that some replicates significant at 2 alpha have the
# wrong sign and the identity is tested in both directions.
small_two_group <- function() build_spec(list(
  name = "dir", seed = 5, design_kind = "between", n_subject = 20,
  factor_name = "group", lev1 = "a", lev2 = "b", intercept = 0, effect = 0.2,
  family = "gaussian", resp_name = "y", sigma = 1))

test_that("one-sided power_design() at alpha is the two-sided test at 2 alpha in one direction", {
  spec <- small_two_group()
  two <- power_design(spec, n_sims = 400, alpha = 0.1)
  greater <- power_design(spec, n_sims = 400, alpha = 0.05, alternative = "greater")
  less <- power_design(spec, n_sims = 400, alpha = 0.05, alternative = "less")
  n_wrong <- round(two$type_s * two$n_significant)
  expect_gt(n_wrong, 0)
  expect_identical(greater$n_significant, two$n_significant - as.integer(n_wrong))
  expect_identical(less$n_significant, as.integer(n_wrong))
  expect_equal(greater$power, two$power * (1 - two$type_s))
  # Type S and Type M are taken over the replicates significant in the tested direction. Against
  # a positive effect, every one of them has the right sign under "greater" and the wrong one
  # under "less", and the exaggeration ratios of the two directions make up the two-sided set.
  expect_identical(greater$type_s, 0)
  expect_identical(less$type_s, 1)
  expect_equal(greater$type_m * greater$n_significant + less$type_m * less$n_significant,
               two$type_m * two$n_significant)
  # The estimates do not depend on the test.
  expect_identical(greater$mean_estimate, two$mean_estimate)
})

test_that("power_design() records its alternative, and two-sided stays the default", {
  spec <- small_two_group()
  default <- power_design(spec, n_sims = 30)
  expect_identical(default$alternative, "two.sided")
  expect_identical(default, power_design(spec, n_sims = 30, alternative = "two.sided"))
  expect_identical(power_design(spec, n_sims = 30, alternative = "less")$alternative, "less")
  expect_error(power_design(spec, n_sims = 2, alternative = "both"))
})

test_that("each two-group replicate's one-sided p is half its two-sided p in its direction", {
  spec <- small_two_group()
  seeds <- replicate_seeds(spec$seed, 6)
  at <- function(i, alternative)
    pilotr:::.power_design_rep(i, spec = spec, seeds = seeds, yname = "y", fname = "group",
                               lev0 = "a", lev1 = "b", alternative = alternative)
  for (i in seq_along(seeds)) {
    two <- at(i, "two.sided")
    toward <- if (two[1] > 0) "greater" else "less"
    away <- setdiff(c("greater", "less"), toward)
    expect_identical(at(i, toward)[2], two[2] / 2, info = i)
    expect_equal(at(i, away)[2], 1 - two[2] / 2, info = i)
  }
})

# The p-values of one replicate of `spec` under each alternative, from the replicate loop that
# power_mixed() runs, and from the fitter the response family chooses.
replicate_ps <- function(spec, i = 1L) {
  seeds <- replicate_seeds(spec$seed, i)
  at <- function(alternative)
    pilotr:::.design_rep(i, spec, seeds, prep = pilotr:::.default_prep(spec),
                         formula = model_formula(spec), fnames = pilotr:::.default_focal(spec),
                         test = TRUE, auto = TRUE, alternative = alternative)
  list(two = at("two.sided"), greater = at("greater"), less = at("less"))
}

crossed_gaussian <- function(effect = 0.3, seed = 1) build_spec(list(
  name = "dirmix", seed = seed, design_kind = "within", include_items = TRUE,
  n_subject = 10, n_item = 8, factor_name = "cond", lev1 = "a", lev2 = "b",
  intercept = 6, effect = effect, subj_int_sd = 0.12, subj_slope_sd = 0,
  item_int_sd = 0.08, item_slope_sd = 0, family = "gaussian", resp_name = "", sigma = 0.3))

crossed_bernoulli <- function() build_spec(list(
  name = "dirglmm", seed = 1, design_kind = "within", include_items = TRUE,
  n_subject = 12, n_item = 10, factor_name = "cond", lev1 = "a", lev2 = "b",
  intercept = 1, effect = 0.5, subj_int_sd = 0.5, subj_slope_sd = 0,
  item_int_sd = 0.3, item_slope_sd = 0, family = "bernoulli", resp_name = ""))

test_that("the one-sided p comes from the same statistic in every fitter", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  # Satterthwaite's t from lmerTest, the t of lm(), and the Wald z of glm() and glmer().
  designs <- list(
    "lmerTest::lmer, Satterthwaite t" = crossed_gaussian(),
    "stats::lm, t" = load_spec(pilotr_example("between_2group_gaussian")),
    "stats::glm (poisson), Wald z" = load_spec(pilotr_example("poisson_counts_between")),
    "lme4::glmer (binomial), Wald z" = crossed_bernoulli())
  for (fitter in names(designs)) {
    r <- replicate_ps(designs[[fitter]])
    tested <- names(r$two$p)[!is.na(r$two$p)]
    expect_gt(length(tested), 0L)
    for (f in tested) {
      info <- paste(fitter, f)
      toward <- if (r$two$est[[f]] > 0) "greater" else "less"
      away <- setdiff(c("greater", "less"), toward)
      expect_identical(r[[toward]]$p[[f]], r$two$p[[f]] / 2, info = info)
      expect_equal(r[[away]]$p[[f]], 1 - r$two$p[[f]] / 2, info = info)
      # The estimate and its standard error do not depend on the test.
      expect_identical(r[[toward]]$est[[f]], r$two$est[[f]], info = info)
      expect_identical(r[[toward]]$se[[f]], r$two$se[[f]], info = info)
    }
  }
})

test_that("one-sided power_mixed() at alpha is the two-sided test at 2 alpha in one direction", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  # No effect is simulated, and the test is set against a positive true value, so that over 16
  # replicates both directions reach significance: 2 in the true direction and 3 in the other.
  spec <- crossed_gaussian(effect = 0, seed = 4)
  focal <- c(effect = 0.05)
  two <- power_mixed(spec, focal = focal, n_sims = 16, alpha = 0.2)
  greater <- power_mixed(spec, focal = focal, n_sims = 16, alpha = 0.1, alternative = "greater")
  less <- power_mixed(spec, focal = focal, n_sims = 16, alpha = 0.1, alternative = "less")
  expect_gt(greater$n_significant[["effect"]], 0L)
  expect_gt(less$n_significant[["effect"]], 0L)
  expect_identical(greater$n_significant + less$n_significant, two$n_significant)
  n_wrong <- round(two$type_s[["effect"]] * two$n_significant[["effect"]])
  expect_identical(less$n_significant[["effect"]], as.integer(n_wrong))
  expect_equal(greater$power, two$power * (1 - two$type_s))
  expect_identical(unname(c(greater$type_s, less$type_s)), c(0, 1))
  # The fit counts and estimates do not depend on the test.
  expect_identical(greater[c("n_returned", "n_converged", "n_singular", "n_warning",
                             "mean_estimate")],
                   two[c("n_returned", "n_converged", "n_singular", "n_warning",
                         "mean_estimate")])
})

test_that("the fitter names a one-sided test, and the default stays two-sided", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  spec <- load_spec(pilotr_example("between_2group_gaussian"))
  expect_identical(power_mixed(spec, n_sims = 2)$fitter, "stats::lm, t")
  greater <- power_mixed(spec, n_sims = 2, alternative = "greater")
  expect_identical(greater$fitter, "stats::lm, t, one-sided (greater)")
  out <- capture.output(print(greater))
  expect_length(grep("one-sided (greater)", out, fixed = TRUE), 1L)
  expect_error(power_mixed(spec, n_sims = 2, alternative = "both"))
  # A sweep passes the alternative to every point, and its rows name the test.
  curve <- power_curve_mixed(spec, subject_ns = c(20, 30), n_sims = 2, alternative = "less")
  expect_identical(curve$fitter, rep("stats::lm, t, one-sided (less)", 2))
  expect_identical(curve$n_significant[2],
                   sweep_spec(spec, "units$subject$n", c(20, 30), power_mixed, n_sims = 2,
                              alternative = "less", .name = "n_subject")$n_significant[2])
  # precision_design() builds an interval, which has no direction, and its fitter is unchanged.
  expect_identical(precision_design(spec, rope = 1, n_sims = 2)$fitter, "stats::lm, Wald z")
})
