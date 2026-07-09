# Parallel execution must reproduce the serial stream exactly. Every replicate seeds the
# shared RNG from its own index, so the worker count can change only the wall-clock time,
# never a number. Cluster tests are skipped on CRAN to respect its two-core limit and to
# avoid PSOCK flakiness on the check farms.

gaussian_between_par <- function(seed = 11) {
  build_spec(list(
    name = "par", seed = seed, n_subject = 24, design_kind = "between",
    factor_name = "group", lev1 = "control", lev2 = "treatment",
    intercept = 100, effect = 5, family = "gaussian", resp_name = "", sigma = 10))
}

crossed_within_par <- function() {
  build_spec(list(
    name = "parmix", seed = 3, n_subject = 12, n_item = 8, design_kind = "within",
    include_items = TRUE, factor_name = "cond", lev1 = "a", lev2 = "b",
    intercept = 6, effect = 0.05, subj_int_sd = 0.12, subj_slope_sd = 0.04, subj_corr = 0.2,
    item_int_sd = 0.08, item_slope_sd = 0.02, item_corr = -0.1,
    family = "shifted_lognormal", resp_name = "", sigma = 0.3, shift = 200))
}

test_that("workers is validated as a single positive whole number", {
  spec <- gaussian_between_par()
  expect_error(power_design(spec, n_sims = 2, workers = 0), "workers")
  expect_error(power_design(spec, n_sims = 2, workers = 1.5), "workers")
  expect_error(power_design(spec, n_sims = 2, workers = c(2, 4)), "workers")
  expect_error(power_design(spec, n_sims = 2, workers = "two"), "workers")
})

test_that("power_design with two workers is identical to a serial run", {
  skip_on_cran()
  spec <- gaussian_between_par()
  serial <- power_design(spec, n_sims = 8, workers = 1)
  parallel2 <- power_design(spec, n_sims = 8, workers = 2)
  expect_identical(serial, parallel2)
})

test_that("power_mixed and power_curve_mixed with two workers are identical to serial runs", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  spec <- crossed_within_par()
  serial <- power_mixed(spec, n_sims = 4, workers = 1)
  parallel2 <- power_mixed(spec, n_sims = 4, workers = 2)
  expect_identical(serial, parallel2)

  curve_serial <- power_curve_mixed(spec, subject_ns = c(12, 14), n_sims = 3, workers = 1)
  curve_parallel <- power_curve_mixed(spec, subject_ns = c(12, 14), n_sims = 3, workers = 2)
  expect_identical(curve_serial, curve_parallel)
})

test_that("precision_design with two workers is identical to a serial run", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  spec <- crossed_within_par()
  serial <- precision_design(spec, focal = c(effect = 0.05), rope = 0.02,
                             n_sims = 4, workers = 1)
  parallel2 <- precision_design(spec, focal = c(effect = 0.05), rope = 0.02,
                                n_sims = 4, workers = 2)
  expect_identical(serial, parallel2)
})
