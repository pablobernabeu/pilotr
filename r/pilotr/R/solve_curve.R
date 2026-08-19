# Solving a simulated design curve for the swept value that meets a target.
#
# The curve functions report a decision rate at each swept value and stop there, leaving the
# reader to judge the crossing from a plot. That judgement is made over points whose Monte Carlo
# intervals overlap, and what it yields is a bare number, the one that goes into a
# preregistration, carrying no record of how wide the range of values compatible with the
# simulation was. The functions here estimate the crossing instead, and report an interval with
# it.
#
# The model is a binomial regression of the rate on the swept value, weighted by the replicate
# count behind each point, with a probit link. The link is not a convenience. A power curve is a
# normal tail probability, and under the normal approximation to a two-group comparison the
# probit of power is exactly linear in the square root of the sample size, so this
# parameterisation estimates two coefficients of a curve the design analysis already implies
# rather than bending a general sigmoid to fit.
#
# The choice was checked rather than assumed, by tools/calibration/solve_curve_calibration.R,
# which produces every figure quoted here and writes them to a file beside itself. Against
# stats::power.t.test(), over twelve combinations of effect size and target power on each of three
# grid shapes at 400 replicates a point, the probit solved to a mean absolute error of 2.75%,
# 2.18% and 3.70%, the logit to 3.27%, 2.54% and 4.28%. The two separate most clearly on the grid
# that sweeps the widest range: reaching a power of 0.995, the logit's Pearson chi-square per
# degree of freedom rose to 1.65 against 1.15 for the probit, and its intervals covered the
# analytic answer nine times in twelve where the probit covered twelve. The margin is not large,
# because at this replicate count most of what separates a solved size from the analytic one is
# Monte Carlo noise in the curve rather than the link fitted to it. Rerun at the seed base 77000,
# which the calibration script takes as its third argument, and the probit still leads on the
# tight grid and by more on the tall one, while the coarse grid swaps: four widely spaced points
# do not tell the two links apart.
#
# The inversion and its interval are the delta method that MASS::dose.p() applies to a fitted
# glm: the solved point on the fitted scale is (link(target) - intercept) / slope, and its
# standard error follows from the gradient of that expression in the two coefficients. The test
# suite checks the agreement with dose.p() where MASS is available.
#
# Fieller's exact interval for a ratio (Fieller, 1954) is the obvious alternative, and it was
# tried. It does not separate from the delta method on this evidence. Over the same 36 checks it
# covered the analytic answer 36 times against the delta method's 35, a difference of one check,
# at a mean width of 9.48 subjects against 9.38, and it was bounded every time. These curves
# determine the slope well enough that the two intervals nearly coincide, so nothing here argues
# for one over the other on coverage. The delta method is kept because it is the interval
# MASS::dose.p() reports, which gives the test suite an independent implementation to check
# against, and because the one condition under which Fieller's would differ, a slope too poorly
# determined for a bounded interval, is the condition the slope refusal below already rejects
# outright. A one-check difference in coverage is noise at this replicate count, and the rerun at
# seed base 77000 above puts the two level at 34 each. Anyone reopening the question should rerun
# the calibration before arguing from it.
#
# One departure from dose.p()'s default: where the two-parameter model does not describe the
# curve, the covariance is scaled by Pearson's chi-square over its degrees of freedom, the
# heterogeneity factor of probit analysis (Finney, 1971). The factor is floored at one, so it
# only ever widens the interval. Without it, a curve the model fits badly still reports the
# narrow interval its replicate counts alone imply, which is the overconfidence this function
# exists to remove.
#
# The regression is written out rather than handed to glm(), because the Python twin has no glm
# to hand it to and the two engines have to agree. With an intercept and one slope the weighted
# normal equations are a two-by-two solve in closed form, small enough to mirror line for line,
# and the reductions are spelled out as loops for the reason core.R gives: base R's sum() and
# CPython's sum() are each more accurate than a plain double fold, and in different ways.
#
# Extrapolation is the one thing these functions must not do. A curve that does not reach the
# target within the range it swept is refused rather than extended, and a fit that solves outside
# that range is refused as well.
#
# Agreement with the Python twin is checked by tools/parity/solve_cross.py, at a relative
# tolerance rather than bit for bit. The fit calls exp() and the normal distribution function at
# every iteration, and IEEE-754 fixes the rounding of neither, so the two languages' maths
# libraries put the last bits in different places even though the arithmetic between them is
# written to be identical.

# Rate columns and replicate-count columns are looked for in this order. The names are those the
# curve functions already use, so a curve goes in unaltered.
.SOLVE_Y <- c("power", "p_meaningful")
.SOLVE_N <- c("n_returned", "n_converged", "n_sims")
.SOLVE_EFFECT <- c("effect", "param")
.SOLVE_TRANSFORMS <- c("sqrt", "identity", "log")
.SOLVE_MAXIT <- 100L
.SOLVE_TOL <- 1e-11
# 1 / sqrt(2 * pi). The normal density is written out with it rather than taken from
# stats::dnorm(), so that the Python twin evaluates the same expression and the only difference
# left between the two is the rounding of exp().
.SOLVE_INV_SQRT_2PI <- 0.3989422804014327

# Every number that reaches a refusal message goes through this, so that the two engines format
# it the same way. Both call the C library's %g on a double.
.solve_num <- function(v) sprintf("%g", v)

# The scale the fit is linear on, and its inverse. The inverse is applied to the solved point and
# to both interval bounds, so it has to be monotone across everything the fit can return; under
# "sqrt" that means flooring the fitted scale at zero before squaring, since a negative square
# root is not a swept value the curve ever visited.
.solve_forward <- function(x, transform)
  switch(transform, sqrt = sqrt(x), identity = x, log = log(x))

.solve_back <- function(u, transform)
  switch(transform, sqrt = { u <- max(u, 0); u * u }, identity = u, log = exp(u))

# Weighted least squares of z on (1, u), returned with the unscaled covariance of the two
# coefficients. The binomial dispersion is one, so this covariance is what the delta method
# needs before any heterogeneity factor is applied to it.
.solve_wls <- function(u, z, w) {
  s0 <- 0; s1 <- 0; s2 <- 0; t0 <- 0; t1 <- 0
  for (i in seq_along(u)) {
    wi <- w[i]
    s0 <- s0 + wi
    s1 <- s1 + wi * u[i]
    s2 <- s2 + wi * u[i] * u[i]
    t0 <- t0 + wi * z[i]
    t1 <- t1 + wi * u[i] * z[i]
  }
  det <- s0 * s2 - s1 * s1
  if (!is.finite(det) || det <= 0) return(NULL)
  list(b0 = (s2 * t0 - s1 * t1) / det, b1 = (s0 * t1 - s1 * t0) / det,
       v00 = s2 / det, v01 = -s1 / det, v11 = s0 / det)
}

# Iteratively reweighted least squares for the binomial probit model with prior weights `m`, the
# number of replicates behind each rate. The starting rates are shrunk towards a half by half a
# replicate, which is what glm() does, so that a rate of exactly zero or one does not start the
# iteration at an infinity.
.solve_irls <- function(u, y, m) {
  p0 <- (m * y + 0.5) / (m + 1)
  eta <- vapply(p0, as241, numeric(1))
  fit <- NULL; prev <- NULL
  for (it in seq_len(.SOLVE_MAXIT)) {
    mu <- stats::pnorm(eta)
    dmu <- exp(-0.5 * eta * eta) * .SOLVE_INV_SQRT_2PI
    v <- mu * (1 - mu)
    # A fitted rate pinned at zero or one has neither a variance to divide by nor a derivative,
    # so the point carries no weight at that step rather than contributing an infinity.
    w <- ifelse(v > 0 & dmu > 0, m * dmu * dmu / v, 0)
    z <- ifelse(dmu > 0, eta + (y - mu) / dmu, eta)
    fit <- .solve_wls(u, z, w)
    if (is.null(fit)) return(NULL)
    b <- c(fit$b0, fit$b1)
    eta <- fit$b0 + fit$b1 * u
    if (!is.null(prev) && max(abs(b - prev)) <= .SOLVE_TOL * (1 + max(abs(b)))) {
      fit$dispersion <- .solve_dispersion(eta, y, m)
      return(fit)
    }
    prev <- b
  }
  NULL
}

# Pearson's chi-square over the residual degrees of freedom, floored at one: the heterogeneity
# factor of probit analysis (Finney, 1971). A curve the two-parameter model does not describe
# leaves residuals larger than the replicate counts alone account for, and scaling the covariance
# by this factor widens the interval to say so. It never narrows one, because a factor below one
# would claim more precision than the replicates behind the curve support.
.solve_dispersion <- function(eta, y, m) {
  df <- length(y) - 2L
  if (df < 1L) return(1)
  mu <- stats::pnorm(eta)
  chi <- 0
  for (i in seq_along(y)) {
    v <- mu[i] * (1 - mu[i])
    if (v > 0) chi <- chi + m[i] * (y[i] - mu[i]) * (y[i] - mu[i]) / v
  }
  max(chi / df, 1)
}

# Resolve a column named by the caller, or the first of `auto` the curve carries.
.solve_column <- function(curve, given, auto, what) {
  cols <- names(curve)
  if (!is.null(given)) {
    if (!(given %in% cols))
      stop(sprintf("`curve` has no column named '%s'. Its columns are: %s.",
                   given, paste(cols, collapse = ", ")), call. = FALSE)
    return(given)
  }
  hit <- auto[auto %in% cols]
  if (!length(hit))
    stop(sprintf("`curve` has no %s column. Name one with `%s`; the columns recognised automatically are %s.",
                 what$label, what$arg, paste(sprintf("'%s'", auto), collapse = ", ")), call. = FALSE)
  hit[1]
}

#' Solve a simulated design curve for the value that meets a target
#'
#' Take the curve a sweep has already produced, fit the decision rate against the swept value,
#' and solve for the value at which the rate meets a target. The solved value comes with a
#' confidence interval, because a point read off a simulated curve without one repeats the
#' overconfidence that design analysis exists to expose.
#'
#' @details
#' The input is the data frame [power_curve_mixed()], [precision_curve()] or [sweep_spec()]
#' returns, used as it stands. The swept value is taken from the leading column, which is where
#' all three put it, and the rate from `power` or `p_meaningful`, whichever the curve carries.
#' Each rate is a proportion over a known number of replicates, and that count, read from
#' `n_returned`, `n_converged` or `n_sims`, weights the fit: a rate over 200 replicates should
#' count for more than a rate over 20.
#'
#' The fit is a binomial regression with a probit link, and the solved value is the swept value
#' at which the fitted rate equals `target`. The probit is chosen because power is a normal tail
#' probability: under the normal approximation to a two-group comparison, the probit of power is
#' linear in the square root of the sample size, so the model has the shape a design analysis
#' already implies. Measured against [stats::power.t.test()] across twelve combinations of effect
#' size and target power on each of three grid shapes, at 400 replicates a point, the solved
#' sample size fell within 2.9% of the analytic answer on average, against 3.4% for a logit
#' fitted the same way.
#'
#' The interval is the delta-method interval of `MASS::dose.p()`, computed on the scale named by
#' `transform` and mapped back, so it is symmetric on that scale rather than on the natural one.
#' That asymmetry is the honest shape: at the top of a power curve a given change in rate costs
#' far more sample size than the same change lower down. Where the two-parameter model does not
#' describe the curve, the interval is widened by the heterogeneity factor of probit analysis,
#' Pearson's chi-square over its degrees of freedom (Finney, 1971), reported as `dispersion`. It
#' is floored at 1, so a well-fitting curve is left alone and a badly-fitting one cannot report a
#' narrower interval than its own residuals justify.
#'
#' The default `transform` of `"sqrt"` suits a sample-size axis, where a rate rises with the
#' square root of the sample size rather than with the sample size itself. Sweep something else,
#' an effect size or a random-effect standard deviation, and `"identity"` is usually right.
#'
#' Nothing here extrapolates. A curve whose rates do not straddle the target is refused, with the
#' range it did cover reported, and so is a fit that solves outside the swept range. A curve
#' whose fitted slope cannot be told from zero is refused too: the crossing is then compatible
#' with any value at all, and an interval that said otherwise would be false.
#'
#' What the interval covers is the Monte Carlo uncertainty of the fit, not the gap between the
#' fitted shape and the true curve. Across 36 checks against [stats::power.t.test()] at 400
#' replicates a point, the solved size sat within 2.9% of the analytic answer on average and
#' within 6.9% at worst, and the nominal 95% interval covered the analytic value 35 times out of
#' 36. Nearly all of that error is the Monte Carlo noise the interval is describing, and it falls
#' with the square root of the replicate count. Raise the count far enough and the interval
#' narrows onto a fitted shape that is still slightly the wrong shape, so replicates alone do not
#' make a solved size arbitrarily accurate. The remedies are a finer grid, more replicates, or a
#' design with a closed form to check against.
#'
#' @param curve A curve, as returned by [power_curve_mixed()], [precision_curve()] or
#'   [sweep_spec()]: one row per swept value, with the swept value, a decision rate, and the
#'   number of replicates behind it.
#' @param target The decision rate to solve for, strictly between 0 and 1.
#' @param x Name of the column holding the swept value. `NULL`, the default, takes the leading
#'   column.
#' @param y Name of the column holding the decision rate. `NULL`, the default, takes `power` or
#'   `p_meaningful`, whichever is present.
#' @param n The number of replicates behind each rate, either the name of a column or a numeric
#'   value. `NULL`, the default, takes `n_returned`, `n_converged` or `n_sims`, whichever is
#'   present.
#' @param effect Which focal effect to solve for, when the curve holds more than one. Matched
#'   against the `effect` or `param` column. `NULL`, the default, uses every row, which is
#'   correct only when the curve holds one effect.
#' @param transform The scale the swept value is fitted on: `"sqrt"` (the default, for a sample
#'   size), `"identity"` or `"log"`.
#' @param level Confidence level for the reported interval.
#' @return A list with elements `value` (the solved swept value), `lo` and `hi` (its confidence
#'   bounds), `level`, `target`, `se` (the delta-method standard error on the fitted scale, the
#'   scale on which the interval is symmetric), `dispersion` (the heterogeneity factor applied,
#'   1 where the model fits), `x` and `y` (the columns used), `transform`, `intercept` and `slope`
#'   (the fitted coefficients), `n_points` (the number of curve points the fit used), and `x_min`
#'   and `x_max` (the swept range). A bound falling outside that range is not an error but a
#'   message: the sweep was too narrow to pin the value down, and should be widened. A
#'   `dispersion` well above 1 says the curve is not the shape the model assumes, so the solve
#'   deserves a wider grid or more replicates rather than trust.
#' @references Fieller, E. C. (1954). Some problems in interval estimation. \emph{Journal of the
#'   Royal Statistical Society: Series B}, 16(2), 175-185.
#'   \doi{10.1111/j.2517-6161.1954.tb00159.x}
#'
#'   Finney, D. J. (1971). \emph{Probit analysis} (3rd ed.). Cambridge University Press.
#' @examples
#' spec <- build_spec(list(name = "s", seed = 1, design_kind = "between", n_subject = 40,
#'   factor_name = "group", lev1 = "a", lev2 = "b", intercept = 0, effect = 0.7,
#'   family = "gaussian", resp_name = "score", sigma = 1))
#' # n_sims is small so the example runs quickly. Use 200 or more for real planning.
#' curve <- sweep_spec(spec, "units$subject$n", c(20, 40, 60, 80), power_design, n_sims = 50)
#' solved <- solve_curve(curve, target = 0.8)
#' unlist(solved[c("value", "lo", "hi")])
#'
#' # This design has an analytic answer to check against, in total subjects across the two
#' # groups. It falls inside the interval, which fifty replicates a point make a wide one.
#' 2 * stats::power.t.test(delta = 0.7, sd = 1, power = 0.8)$n
#' @seealso [target_n()] for the sample-size case, and [power_curve_mixed()],
#'   [precision_curve()] and [sweep_spec()] for the curves this consumes.
#' @export
solve_curve <- function(curve, target, x = NULL, y = NULL, n = NULL, effect = NULL,
                        transform = "sqrt", level = 0.95) {
  if (!is.data.frame(curve) || !ncol(curve) || !nrow(curve))
    stop("`curve` must be a table of curve points, with one row per swept value.", call. = FALSE)
  .solve_probability(target, "target")
  .solve_probability(level, "level")
  if (!is.character(transform) || length(transform) != 1L ||
      !(transform %in% .SOLVE_TRANSFORMS))
    stop(sprintf("`transform` must be one of %s.",
                 paste(sprintf("'%s'", .SOLVE_TRANSFORMS), collapse = ", ")), call. = FALSE)

  xname <- if (is.null(x)) names(curve)[1] else
    .solve_column(curve, x, character(0), list(label = "swept-value", arg = "x"))
  yname <- .solve_column(curve, y, .SOLVE_Y, list(label = "decision-rate", arg = "y"))

  keep <- .solve_rows(curve, effect)
  xv <- .solve_numeric(curve[[xname]][keep], xname, "swept-value")
  yv <- .solve_numeric(curve[[yname]][keep], yname, "decision-rate")
  mv <- .solve_weights(curve, n, keep)

  ok <- is.finite(xv) & is.finite(yv) & is.finite(mv) & mv > 0
  xv <- xv[ok]; yv <- yv[ok]; mv <- mv[ok]
  k <- length(xv)
  if (k < 3L)
    stop(sprintf("solving a curve needs at least 3 swept values with a rate and a replicate count; this curve has %d.",
                 k), call. = FALSE)
  if (anyDuplicated(xv))
    stop("the curve has more than one row at the same swept value. Select a single focal effect with `effect`, or subset the curve before solving.",
         call. = FALSE)
  bad <- yv[yv < 0 | yv > 1]
  if (length(bad))
    stop(sprintf("the column '%s' holds a value of %s, which is not a probability; solve_curve() inverts a rate between 0 and 1.",
                 yname, .solve_num(bad[1])), call. = FALSE)
  if (max(yv) == min(yv))
    stop(sprintf("the rate is %s at every swept value, so the curve has no trend to invert.",
                 .solve_num(yv[1])), call. = FALSE)
  if (transform == "sqrt" && min(xv) < 0)
    stop(sprintf("the 'sqrt' transform is not defined for a swept value of %s. Use transform = \"identity\" for an axis that is not a sample size.",
                 .solve_num(min(xv))), call. = FALSE)
  if (transform == "log" && min(xv) <= 0)
    stop(sprintf("the 'log' transform is not defined for a swept value of %s. Use transform = \"identity\" for an axis that is not a sample size.",
                 .solve_num(min(xv))), call. = FALSE)
  if (target < min(yv) || target > max(yv))
    stop(sprintf("the curve does not reach a %s of %s within its swept range. Over swept values from %s to %s the rate runs from %s to %s, and solving would extrapolate beyond the values simulated.",
                 yname, .solve_num(target), .solve_num(min(xv)), .solve_num(max(xv)),
                 .solve_num(min(yv)), .solve_num(max(yv))), call. = FALSE)

  u <- .solve_forward(xv, transform)
  fit <- .solve_irls(u, yv, mv)
  if (is.null(fit))
    stop(sprintf("the weighted fit to this curve did not settle within %d iterations, so the curve cannot be inverted.",
                 .SOLVE_MAXIT), call. = FALSE)

  z <- as241((1 + level) / 2)
  se_slope <- sqrt(fit$dispersion * fit$v11)
  if (!is.finite(se_slope) || abs(fit$b1) <= z * se_slope)
    stop(sprintf("the fitted slope (%s, standard error %s) cannot be told from zero at the %s level, so the curve does not determine a value; every swept value is compatible with the target.",
                 .solve_num(fit$b1), .solve_num(se_slope), .solve_num(level)), call. = FALSE)

  et <- as241(target)
  us <- (et - fit$b0) / fit$b1
  g0 <- -1 / fit$b1
  g1 <- -us / fit$b1
  se <- sqrt(fit$dispersion *
               (g0 * g0 * fit$v00 + 2 * g0 * g1 * fit$v01 + g1 * g1 * fit$v11))
  value <- .solve_back(us, transform)
  if (value < min(xv) || value > max(xv))
    stop(sprintf("the solved value %s falls outside the swept range %s to %s, so reporting it would extrapolate beyond the values simulated.",
                 .solve_num(value), .solve_num(min(xv)), .solve_num(max(xv))), call. = FALSE)

  list(value = value, lo = .solve_back(us - z * se, transform),
       hi = .solve_back(us + z * se, transform),
       level = level, target = target, se = se, dispersion = fit$dispersion,
       x = xname, y = yname, transform = transform,
       intercept = fit$b0, slope = fit$b1, n_points = k,
       x_min = min(xv), x_max = max(xv))
}

#' Solve a power curve for the sample size that reaches a target power
#'
#' The sample-size case of [solve_curve()], and the number a power analysis is usually run to
#' obtain. Takes the curve a sweep over sample size has produced and returns the size at which
#' power reaches `target`, rounded up to a whole number of units alongside the exact solution.
#'
#' @details
#' Everything [solve_curve()] does applies here, including its refusals: a curve that never
#' reaches the target within the sizes it swept is refused rather than extrapolated, and the
#' reported interval can extend past the largest size simulated, which means the sweep was too
#' narrow to settle the question.
#'
#' The whole-number fields round up rather than to nearest, because a design cannot recruit a
#' fraction of a subject and rounding down would leave the study short of the target it was
#' sized for.
#'
#' @param curve A power curve, as returned by [power_curve_mixed()] or by [sweep_spec()] over
#'   `units$subject$n`.
#' @param target The power to reach. Defaults to 0.8, the convention this package's plots draw a
#'   line at.
#' @param ... Further arguments passed to [solve_curve()], such as `effect` to pick one focal
#'   effect out of a curve holding several, or `level` for the interval.
#' @return The list [solve_curve()] returns, with `n`, `n_lo` and `n_hi` added: `value`, `lo` and
#'   `hi` rounded up to whole numbers.
#' @examples
#' spec <- build_spec(list(name = "s", seed = 1, design_kind = "between", n_subject = 40,
#'   factor_name = "group", lev1 = "a", lev2 = "b", intercept = 0, effect = 0.7,
#'   family = "gaussian", resp_name = "score", sigma = 1))
#' # n_sims is small so the example runs quickly. Use 200 or more for real planning.
#' curve <- sweep_spec(spec, "units$subject$n", c(20, 40, 60, 80), power_design, n_sims = 50)
#' solved <- target_n(curve)
#' unlist(solved[c("n", "n_lo", "n_hi")])
#' @seealso [solve_curve()], which this wraps, and [power_curve_mixed()] for the curve.
#' @export
target_n <- function(curve, target = 0.8, ...) {
  out <- solve_curve(curve, target = target, ...)
  c(out, list(n = ceiling(out$value), n_lo = ceiling(out$lo), n_hi = ceiling(out$hi)))
}

# A single probability strictly inside the unit interval. The two checks are separate so that a
# missing or non-numeric argument is never formatted into the message, which R would render as
# "NA" and Python as "nan".
.solve_probability <- function(v, arg) {
  if (!is.numeric(v) || length(v) != 1L || !is.finite(v))
    stop(sprintf("`%s` must be a single number strictly between 0 and 1.", arg), call. = FALSE)
  if (v <= 0 || v >= 1)
    stop(sprintf("`%s` must be strictly between 0 and 1; got %s.", arg, .solve_num(v)),
         call. = FALSE)
}

# The rows to fit, after selecting a focal effect if one was named.
.solve_rows <- function(curve, effect) {
  if (is.null(effect)) return(seq_len(nrow(curve)))
  col <- .SOLVE_EFFECT[.SOLVE_EFFECT %in% names(curve)]
  if (!length(col))
    stop("`curve` has no 'effect' or 'param' column to select a focal effect from.", call. = FALSE)
  have <- as.character(curve[[col[1]]])
  keep <- which(have == as.character(effect))
  if (!length(keep))
    stop(sprintf("`curve` holds no focal effect named '%s'. It holds: %s.",
                 effect, paste(unique(have), collapse = ", ")), call. = FALSE)
  keep
}

.solve_numeric <- function(v, nm, what) {
  if (!is.numeric(v))
    stop(sprintf("the %s column '%s' is not numeric.", what, nm), call. = FALSE)
  as.numeric(v)
}

# The replicate count behind each rate: a column name, a numeric value recycled over the rows, or
# the first count column the curve carries.
.solve_weights <- function(curve, n, keep) {
  if (is.numeric(n)) {
    if (!length(n)) stop("`n` must hold at least one replicate count.", call. = FALSE)
    return(rep(as.numeric(n), length.out = length(keep)))
  }
  nm <- .solve_column(curve, n, .SOLVE_N, list(label = "replicate-count", arg = "n"))
  .solve_numeric(curve[[nm]][keep], nm, "replicate-count")
}
