# Tests for solving a simulated design curve. The curves are written out rather than simulated,
# so that a change to the generative core cannot move a number here, and so that the Python
# suite can hold exactly the same tables in test_solve.py.
#
# analytic_curve() is the power of a two-sample t-test at a standardised effect of 0.5, taken
# from stats::power.t.test() at each total sample size and rounded to four places. It is the
# external reference the solver is checked against: the same function reports that 0.80 power
# needs 127.5315 subjects in total and 0.90 needs 170.0626, and the solve has to land near those
# without ever having been told them.

analytic_curve <- function(n_sims = 2000) {
  data.frame(n_subject = c(40, 60, 80, 100, 120, 140, 160, 180),
             power = c(0.3377, 0.4778, 0.5981, 0.6969, 0.7753, 0.8358, 0.8816, 0.9156),
             n_sims = n_sims)
}

ANALYTIC_N80 <- 127.5315
ANALYTIC_N90 <- 170.0626

test_that("the solved sample size matches the analytic answer for a design that has one", {
  cv <- analytic_curve()
  s80 <- solve_curve(cv, target = 0.8)
  s90 <- solve_curve(cv, target = 0.9)
  expect_lt(abs(s80$value - ANALYTIC_N80) / ANALYTIC_N80, 0.02)
  expect_lt(abs(s90$value - ANALYTIC_N90) / ANALYTIC_N90, 0.02)
  # The interval is the point of the exercise, so it has to cover the value it is estimating.
  expect_true(s80$lo <= ANALYTIC_N80 && ANALYTIC_N80 <= s80$hi)
  expect_true(s90$lo <= ANALYTIC_N90 && ANALYTIC_N90 <= s90$hi)
})

test_that("the returned solve describes itself", {
  s <- solve_curve(analytic_curve(), target = 0.8)
  expect_setequal(names(s), c("value", "lo", "hi", "level", "target", "se", "dispersion",
                              "x", "y", "transform", "intercept", "slope", "n_points",
                              "x_min", "x_max"))
  expect_identical(s$x, "n_subject")
  expect_identical(s$y, "power")
  expect_identical(s$transform, "sqrt")
  expect_identical(s$n_points, 8L)
  expect_equal(s$x_min, 40)
  expect_equal(s$x_max, 180)
  expect_lt(s$lo, s$value)
  expect_lt(s$value, s$hi)
  expect_gt(s$slope, 0)
})

test_that("a wider confidence level gives a wider interval", {
  cv <- analytic_curve()
  narrow <- solve_curve(cv, target = 0.8, level = 0.8)
  wide <- solve_curve(cv, target = 0.8, level = 0.99)
  expect_equal(narrow$value, wide$value)
  expect_lt(wide$lo, narrow$lo)
  expect_gt(wide$hi, narrow$hi)
})

test_that("more replicates behind the same rates narrow the interval but leave the point alone", {
  few <- solve_curve(analytic_curve(n_sims = 50), target = 0.8)
  many <- solve_curve(analytic_curve(n_sims = 5000), target = 0.8)
  # A weight common to every point cancels out of the fitted line, though not out of the last
  # bits of it, so the two solves agree to rounding rather than exactly.
  expect_equal(few$value, many$value, tolerance = 1e-12)
  expect_lt(many$hi - many$lo, few$hi - few$lo)
})

test_that("target_n rounds up, because a design cannot recruit a fraction of a subject", {
  cv <- analytic_curve()
  s <- target_n(cv)
  expect_equal(s$target, 0.8)
  expect_identical(s$n, ceiling(s$value))
  expect_identical(s$n_lo, ceiling(s$lo))
  expect_identical(s$n_hi, ceiling(s$hi))
  expect_gte(s$n, s$value)
  expect_lt(s$n - s$value, 1)
})

test_that("the fit agrees with a binomial glm and the inversion with MASS::dose.p", {
  skip_if_not_installed("MASS")
  # glm() takes counts rather than a rate and a weight, so the replicate count is set to make
  # every rate a whole number of successes. Otherwise the two fits would see different data
  # and the comparison would measure the rounding rather than the arithmetic.
  cv <- analytic_curve(n_sims = 10000)
  u <- sqrt(cv$n_subject)
  successes <- round(cv$power * cv$n_sims)
  failures <- cv$n_sims - successes
  fit <- stats::glm(cbind(successes, failures) ~ u, family = stats::binomial(link = "probit"))
  s <- solve_curve(cv, target = 0.8)
  expect_equal(unname(stats::coef(fit)), c(s$intercept, s$slope), tolerance = 1e-10)
  # dose.p() reports the solve on the fitted scale, which here is the square root of the
  # sample size, together with the delta-method standard error solve_curve() reports as `se`.
  dp <- MASS::dose.p(fit, p = 0.8)
  expect_equal(as.numeric(dp), sqrt(s$value), tolerance = 1e-10)
  # The standard error is looser because it is read off the weights of the last iteration, and
  # glm() stops on a relative change in deviance of 1e-8 while solve_curve() runs the
  # coefficients down to 1e-11.
  expect_equal(as.numeric(attr(dp, "SE")), s$se, tolerance = 1e-6)
})

test_that("a precision curve is solved through its own rate column", {
  cv <- data.frame(n_subject = c(15, 30, 60, 100, 140, 180, 220, 260),
                   p_meaningful = c(0.06, 0.14, 0.33, 0.55, 0.71, 0.83, 0.90, 0.94),
                   n_returned = 200)
  s <- solve_curve(cv, target = 0.9)
  expect_identical(s$y, "p_meaningful")
  expect_gt(s$value, 180)
  expect_lt(s$value, 260)
})

test_that("a non-sample-size axis is solved on the identity scale", {
  cv <- data.frame(effect_size = c(0.1, 0.2, 0.3, 0.4, 0.5),
                   power = c(0.11, 0.31, 0.58, 0.79, 0.92), n_sims = 400)
  s <- solve_curve(cv, target = 0.8, transform = "identity")
  expect_identical(s$x, "effect_size")
  expect_gt(s$value, 0.3)
  expect_lt(s$value, 0.5)
})

test_that("a focal effect is selected out of a curve holding several", {
  cv <- data.frame(n_subject = rep(c(20, 40, 60, 80, 100), each = 2),
                   effect = rep(c("a", "b"), 5),
                   power = c(0.20, 0.10, 0.40, 0.18, 0.62, 0.27, 0.78, 0.35, 0.88, 0.44),
                   n_sims = 300)
  s <- solve_curve(cv, target = 0.7, effect = "a")
  expect_equal(s$n_points, 5L)
  expect_gt(s$value, 40)
  expect_lt(s$value, 100)
  expect_error(solve_curve(cv, target = 0.7),
               "the curve has more than one row at the same swept value", fixed = TRUE)
  expect_error(solve_curve(cv, target = 0.7, effect = "z"),
               "`curve` holds no focal effect named 'z'. It holds: a, b.", fixed = TRUE)
})

test_that("the replicate count may be given rather than read off a column", {
  cv <- data.frame(n_subject = c(40, 80, 120, 160), power = c(0.34, 0.60, 0.78, 0.88))
  expect_error(solve_curve(cv, target = 0.7),
               "`curve` has no replicate-count column. Name one with `n`; the columns recognised automatically are 'n_returned', 'n_converged', 'n_sims'.",
               fixed = TRUE)
  s <- solve_curve(cv, target = 0.7, n = 500)
  expect_gt(s$value, 40)
  expect_lt(s$value, 160)
})

# --- refusals ---------------------------------------------------------------------------

test_that("a target outside the unit interval is refused", {
  cv <- analytic_curve()
  for (bad in c(0, 1, 1.5, -0.2)) {
    expect_error(solve_curve(cv, target = bad),
                 sprintf("`target` must be strictly between 0 and 1; got %s.", format(bad)),
                 fixed = TRUE)
  }
  expect_error(solve_curve(cv, target = NA_real_),
               "`target` must be a single number strictly between 0 and 1.", fixed = TRUE)
  expect_error(solve_curve(cv, target = c(0.5, 0.8)),
               "`target` must be a single number strictly between 0 and 1.", fixed = TRUE)
  expect_error(solve_curve(cv, target = 0.8, level = 1),
               "`level` must be strictly between 0 and 1; got 1.", fixed = TRUE)
})

test_that("a curve that never reaches the target is refused, with the range it did cover", {
  expect_error(
    solve_curve(analytic_curve(), target = 0.95),
    "the curve does not reach a power of 0.95 within its swept range. Over swept values from 40 to 180 the rate runs from 0.3377 to 0.9156, and solving would extrapolate beyond the values simulated.",
    fixed = TRUE)
})

test_that("a solve landing outside the swept range is refused rather than extrapolated", {
  # The target is exactly the lowest rate simulated, so the bracket test passes; the fitted
  # curve sits a shade above that point, so the solve falls just short of the smallest size
  # swept, which is outside the range whether it misses by one subject or a hundred.
  expect_error(solve_curve(analytic_curve(), target = 0.3377),
               "falls outside the swept range 40 to 180, so reporting it would extrapolate beyond the values simulated.",
               fixed = TRUE)
})

test_that("a target exactly on the highest rate simulated is solved rather than refused", {
  s <- solve_curve(analytic_curve(), target = 0.9156)
  expect_gte(s$value, 40)
  expect_lte(s$value, 180)
})

test_that("a curve the model does not describe reports a dispersion above 1 and a wider interval", {
  cv <- analytic_curve()
  clean <- solve_curve(cv, target = 0.8)
  expect_equal(clean$dispersion, 1)
  # The same curve with its middle points pushed about, which no two-parameter model can pass
  # through. The point estimate barely moves; the interval has to say that it is less certain.
  rough <- cv
  rough$power <- c(0.3377, 0.5600, 0.5400, 0.7600, 0.7100, 0.8700, 0.8500, 0.9156)
  ragged <- solve_curve(rough, target = 0.8)
  expect_gt(ragged$dispersion, 1)
  expect_gt(ragged$hi - ragged$lo, 3 * (clean$hi - clean$lo))
})

test_that("a curve with too few points is refused", {
  cv <- analytic_curve()
  expect_error(solve_curve(cv[1:2, ], target = 0.5),
               "solving a curve needs at least 3 swept values with a rate and a replicate count; this curve has 2.",
               fixed = TRUE)
  expect_error(solve_curve(cv[1, , drop = FALSE], target = 0.5),
               "solving a curve needs at least 3 swept values with a rate and a replicate count; this curve has 1.",
               fixed = TRUE)
  expect_error(solve_curve(cv[0, ], target = 0.5),
               "`curve` must be a table of curve points, with one row per swept value.",
               fixed = TRUE)
  # A point with no converged fits has no rate, and is dropped before the count is taken.
  thin <- cv[1:4, ]; thin$power[2:3] <- NA_real_
  expect_error(solve_curve(thin, target = 0.5),
               "solving a curve needs at least 3 swept values with a rate and a replicate count; this curve has 2.",
               fixed = TRUE)
})

test_that("a curve with no trend in the rate is refused", {
  flat <- data.frame(n_subject = c(10, 20, 30), power = 0.5, n_sims = 100)
  expect_error(solve_curve(flat, target = 0.5),
               "the rate is 0.5 at every swept value, so the curve has no trend to invert.",
               fixed = TRUE)
})

test_that("a curve whose slope cannot be told from zero is refused", {
  noise <- data.frame(n_subject = c(10, 20, 30, 40, 50),
                      power = c(0.50, 0.52, 0.49, 0.51, 0.50), n_sims = 30)
  expect_error(solve_curve(noise, target = 0.5),
               "cannot be told from zero at the 0.95 level, so the curve does not determine a value; every swept value is compatible with the target.",
               fixed = TRUE)
})

test_that("a transform undefined on the swept values is refused", {
  cv <- data.frame(delta = c(-0.2, 0, 0.2), power = c(0.20, 0.50, 0.85), n_sims = 200)
  expect_error(solve_curve(cv, target = 0.6),
               "the 'sqrt' transform is not defined for a swept value of -0.2. Use transform = \"identity\" for an axis that is not a sample size.",
               fixed = TRUE)
  expect_error(solve_curve(cv, target = 0.6, transform = "log"),
               "the 'log' transform is not defined for a swept value of -0.2. Use transform = \"identity\" for an axis that is not a sample size.",
               fixed = TRUE)
  s <- solve_curve(cv, target = 0.6, transform = "identity")
  expect_gt(s$value, -0.2)
  expect_lt(s$value, 0.2)
})

test_that("a column that is not a rate, or not there at all, is refused", {
  cv <- analytic_curve()
  expect_error(solve_curve(cv, target = 0.8, y = "nope"),
               "`curve` has no column named 'nope'. Its columns are: n_subject, power, n_sims.",
               fixed = TRUE)
  expect_error(solve_curve(data.frame(n = 1:5, width = c(0.9, 0.7, 0.5, 0.4, 0.3)), target = 0.5),
               "`curve` has no decision-rate column. Name one with `y`; the columns recognised automatically are 'power', 'p_meaningful'.",
               fixed = TRUE)
  wide <- data.frame(n = 1:4, power = c(0.1, 0.5, 0.9, 1.4), n_sims = 100)
  expect_error(solve_curve(wide, target = 0.5),
               "the column 'power' holds a value of 1.4, which is not a probability; solve_curve() inverts a rate between 0 and 1.",
               fixed = TRUE)
  expect_error(solve_curve(cv, target = 0.8, transform = "cube"),
               "`transform` must be one of 'sqrt', 'identity', 'log'.", fixed = TRUE)
})

test_that("solve_curve consumes what sweep_spec produces", {
  spec <- build_spec(list(name = "s", seed = 11, design_kind = "between", n_subject = 40,
    factor_name = "group", lev1 = "a", lev2 = "b", intercept = 0, effect = 0.6,
    family = "gaussian", resp_name = "score", sigma = 1))
  cv <- sweep_spec(spec, "units$subject$n", c(20, 40, 60, 80, 100), power_design,
                   n_sims = 200, .name = "n_subject")
  s <- target_n(cv)
  expect_identical(s$x, "n_subject")
  expect_identical(s$y, "power")
  expect_gte(s$n, 20)
  expect_lte(s$n, 100)
})
