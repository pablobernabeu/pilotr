# The grid builder and the variance calibrator are exported and recommended by the sweep
# documentation, but were reached only through the analyses that consume them, so their own
# shapes and their refusals went unchecked. Everything here runs on the variance decomposition
# rather than on a model fit, so no replicate loop is involved.

test_that("design_conditions recycles, prepends the null and carries a base", {
  cond <- design_conditions(cond = c(0.02, 0.05), age = 0.1)
  expect_length(cond, 3L)
  expect_identical(cond[[1]], list(cond = 0, age = 0))   # the shared all-zero condition
  expect_identical(cond[[2]], list(cond = 0.02, age = 0.1))
  expect_identical(cond[[3]], list(cond = 0.05, age = 0.1))

  plain <- design_conditions(effect = c(0.03, 0.06), .null = FALSE)
  expect_length(plain, 2L)
  expect_identical(plain[[1]], list(effect = 0.03))

  based <- design_conditions(effect = 0.03, .null = FALSE, .base = list(age = 0.2))
  expect_identical(based[[1]], list(age = 0.2, effect = 0.03))
})

test_that("design_conditions refuses an unnamed or non-numeric effect", {
  expect_error(design_conditions(c(0.02, 0.05)), "name every effect")
  expect_error(design_conditions(cond = "0.02"), "must be numeric")
})

test_that("sweep_spec walks a grid of coefficient sets by index", {
  spec <- load_spec(pilotr_example("between_2group_gaussian"))
  total <- function(spec, ...) data.frame(total = response_variance(spec)$total)
  out <- sweep_spec(spec, "fixed$coefficients",
                    design_conditions(grp = c(1, 2)), total)
  expect_equal(nrow(out), 3L)
  # The swept values are lists, so the leading column records the grid index instead.
  expect_identical(names(out)[1], "coefficients")
  expect_identical(out$coefficients, 1:3)
  expect_true(all(diff(out$total) > 0))
})

test_that("sweep_spec refuses a path that addresses no field", {
  spec <- load_spec(pilotr_example("between_2group_gaussian"))
  expect_error(sweep_spec(spec, "units$cluster$n", c(10, 20), power_design),
               "does not address a field")
  expect_error(sweep_spec(spec, "units$subject$n", numeric(0), power_design),
               "at least one value")
  expect_error(sweep_spec(spec, "units$subject$n", c(10, 20), "power_design"),
               "must be a function")
})

test_that("calibrate_response reaches the target for a family whose residual moves with eta", {
  # beta is solved numerically, since its residual depends on the linear predictor.
  spec <- load_spec(pilotr_example("beta_proportion"))
  expect_equal(response_variance(calibrate_response(spec, 5, tune = "all"))$total, 5)
  pois <- load_spec(pilotr_example("poisson_counts_between"))
  expect_equal(response_variance(calibrate_response(pois, 20, tune = "all"))$total, 20)
})

test_that("calibrate_response refuses the targets no rescaling can reach", {
  ord <- load_spec(pilotr_example("ordinal_likert_between"))
  expect_error(calibrate_response(ord, 1, tune = "sigma"),
               "has no residual standard deviation to tune")
  # pi^2 / 3 on the latent scale, which no rescaling of the design can move.
  expect_error(calibrate_response(ord, 1, tune = "all"),
               "already meets or exceeds the target")

  gauss <- load_spec(pilotr_example("between_2group_gaussian"))
  expect_error(calibrate_response(gauss, 0.001, tune = "sigma"),
               "no residual standard deviation can reach it")
  expect_error(calibrate_response(gauss, 0, tune = "all"),
               "must be a single positive number")
})

test_that("calibrate_response holds beta fixed when it tunes an exgaussian sigma", {
  spec <- build_spec(list(name = "eg", seed = 1, design_kind = "between",
                          factor_name = "g", lev1 = "a", lev2 = "b", n_subject = 40,
                          intercept = 0, effect = 0.4, family = "exgaussian",
                          resp_name = "", sigma = 0.5, beta = 0.4))
  tuned <- calibrate_response(spec, 1, tune = "sigma")
  expect_equal(tuned$response$beta, 0.4)          # only sigma absorbs the difference
  expect_equal(response_variance(tuned)$total, 1)

  spec$response$beta <- 5
  expect_error(calibrate_response(spec, 1, tune = "sigma"), "leaving nothing for sigma")
})
