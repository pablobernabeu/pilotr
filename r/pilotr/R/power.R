# Simulation-based power and design analysis (Type S and Type M), mirroring power.py.

#' Simulation-based power and design analysis for a two-group Gaussian design
#'
#' Estimate power by repeatedly simulating from the specification and applying a two-sample
#' t-test, alongside the Type S (sign) and Type M (magnitude) design-analysis errors of
#' Gelman and Carlin (2014).
#'
#' @details
#' The t-test treats every row as an independent observation, which is valid only when each
#' subject contributes one row and no rows share a cluster. `power_design()` therefore takes a
#' Gaussian design with exactly one two-level factor between subjects and no item unit, within
#' factor or grouping factor besides `subject`. Continuous predictors and by-subject random
#' effects are allowed, since with one row per subject they vary independently from row to row.
#' Any other design is refused with a message naming what makes its rows correlated.
#' [power_mixed()] fits the model that such a design implies. Take 30 subjects crossed with 20
#' items, with by-subject and residual standard deviations of 1 and a by-item one of 0.3. With
#' no true effect, a t-test of every row is significant in about half of all replicates.
#'
#' A specification with no coefficient for the factor's contrast has a true effect of 0, so
#' `type_s` and `type_m` are `NaN`.
#'
#' Power, Type S and Type M are each estimated from a finite number of replicates, so each is
#' reported with its Monte Carlo error. Power has its standard error and a 95% Wilson interval,
#' one of the intervals Brown, Cai and DasGupta (2001) recommend in place of the Wald interval. It
#' stays between 0 and 1 and does not shrink to a point when power is 0 or 1. Type S and Type M
#' rest on the significant replicates alone, often a handful. Type S has a Wilson interval,
#' because its standard error is zero whenever no sign is wrong. Type M has its standard error,
#' the standard deviation of the exaggeration ratios over the square root of their number.
#'
#' A directional prediction calls for a one-sided test, which `alternative` selects. The estimate
#' is the mean of the second level minus that of the first, whose true value is `true_effect`.
#' `"greater"` tests whether it is positive and `"less"` whether it is negative, by
#' `stats::t.test()` with the same `alternative`. The one-sided p-value is half the two-sided one
#' when the estimate lies in the tested direction. A one-sided test at `alpha` therefore rejects in
#' exactly the replicates in which the two-sided test at `2 * alpha` rejects with the estimate in
#' that direction. Type S and Type M are taken over the replicates significant in the tested
#' direction. When that is the direction of the true effect, Type S is 0 by construction, and
#' Type M is the exaggeration among the replicates that would support the prediction.
#'
#' Results from earlier versions were two-sided. A one-sided power at `alpha`, in the direction of
#' the true effect, can be read off a two-sided result at `2 * alpha`: it is
#' `power * (1 - type_s)`. Type M is recovered only when that result's Type S is 0, since the
#' wrong-signed replicates cannot otherwise be told apart from the others.
#'
#' @param spec A design specification (path or list) for a two-group Gaussian design with one
#'   row per subject (see Details).
#' @param n_sims Number of Monte Carlo replicates. A power estimate carries a Monte Carlo
#'   standard error of about `sqrt(p * (1 - p) / n_sims)`, reported as `power_mcse`, and
#'   `type_s` and `type_m` average over the significant replicates alone, so they settle more
#'   slowly still. At least 200 replicates are advisable for study planning.
#' @param alpha Significance level of the test, two-sided unless `alternative` says otherwise.
#' @param workers Number of local worker processes over which to spread the replicates.
#'   The default of 1 runs serially. Because every replicate takes its own seed from
#'   [replicate_seeds()], any worker count returns results identical to a serial run.
#' @param cl A cluster made by [parallel::makeCluster()] to run the replicates on, in place of one
#'   started for this call, in which case `workers` is ignored. The cluster is left running for its
#'   owner to stop, and through [sweep_spec()] one cluster serves every grid point.
#' @param alternative The test's alternative hypothesis: `"two.sided"`, the default, or
#'   `"greater"` or `"less"` for a one-sided test of a directional prediction (see Details).
#' @return A list with elements `n_sims`, `alpha`, `alternative`, `power`, `power_mcse` (its
#'   Monte Carlo standard error), `power_lo` and `power_hi` (its 95% Wilson interval),
#'   `n_significant`, `true_effect` and `mean_estimate`. Then come `type_s` (sign-error rate
#'   among significant replicates), `type_s_lo` and `type_s_hi` (its 95% Wilson interval),
#'   `type_m` (mean exaggeration ratio among significant replicates) and `type_m_mcse` (its Monte
#'   Carlo standard error). Both design-analysis quantities are `NaN` when no replicate reached
#'   significance and when the true effect is zero, as in the null condition
#'   [design_conditions()] produces: neither is defined without a true value to compare
#'   against, and Type M divides by it. Their interval and standard error are `NaN` then too,
#'   and `type_m_mcse` is also `NaN` when a single replicate reached significance, since a
#'   standard deviation needs two.
#' @references Brown, L. D., Cai, T. T. and DasGupta, A. (2001). Interval estimation for a
#'   binomial proportion. \emph{Statistical Science}, 16(2), 101-133.
#'   \doi{10.1214/ss/1009213286}
#'
#'   Gelman, A. and Carlin, J. (2014). Beyond power calculations: Assessing Type S
#'   (sign) and Type M (magnitude) errors. \emph{Perspectives on Psychological Science},
#'   9(6), 641-651. \doi{10.1177/1745691614551642}
#' @examples
#' spec <- build_spec(list(name = "d", seed = 1, design_kind = "between",
#'   factor_name = "group", lev1 = "a", lev2 = "b", n_subject = 64,
#'   intercept = 100, effect = 5, family = "gaussian", resp_name = "", sigma = 10))
#' # n_sims is small so the example runs quickly. Use 200 or more for real planning.
#' power_design(spec, n_sims = 50)
#' @export
power_design <- function(spec, n_sims = 1000, alpha = 0.05, workers = 1, cl = NULL,
                         alternative = c("two.sided", "greater", "less")) {
  alternative <- match.arg(alternative)
  spec <- .as_spec(spec)
  refusal <- .two_group_refusal(spec)
  if (!is.null(refusal)) stop(refusal, call. = FALSE)
  f <- Filter(function(f) !is.null(f[["between"]]), spec[["factors"]])[[1]]
  fname <- f[["name"]]; lev0 <- f[["levels"]][1]; lev1 <- f[["levels"]][2]
  col <- names(f[["contrasts"]])[1]; vals <- f[["contrasts"]][[col]]
  # A specification may leave the contrast out of `coefficients`, `{}` included, and the
  # simulator then generates no effect, so the true effect is 0. A NULL here used to reach the
  # Type S and Type M guard as a zero-length value and stop the run.
  beta <- spec[["fixed"]][["coefficients"]][[col]]
  if (is.null(beta)) beta <- 0
  true_effect <- beta * (vals[2] - vals[1])
  yname <- spec[["response"]][["name"]]
  seeds <- replicate_seeds(spec[["seed"]], n_sims)

  workers <- .check_workers(workers)
  cl <- .check_cluster(cl)
  if (is.null(cl) && workers > 1L) {
    cl <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }
  res <- .p_lapply(seq_len(n_sims), .power_design_rep, cl = cl, spec = spec, seeds = seeds,
                   yname = yname, fname = fname, lev0 = lev0, lev1 = lev1,
                   alternative = alternative)
  est <- vapply(res, `[[`, numeric(1), 1L)
  pv  <- vapply(res, `[[`, numeric(1), 2L)
  # Under a one-sided test, these are the replicates significant in the tested direction, so Type
  # S and Type M below are taken over them alone.
  sig <- which(pv < alpha)
  n_sig <- length(sig)
  # Type S and Type M are defined relative to a true value, and Type M divides by it, so both
  # stay NaN when the true effect is zero. The alternative was an infinity, alongside a
  # sign-error rate that had quietly become "the estimate is positive". Same rule as
  # power_mixed(), and the null condition design_conditions() recommends is exactly this case.
  usable <- n_sig > 0L && !is.na(true_effect) && true_effect != 0
  wrong <- (est[sig] > 0) != (true_effect > 0)
  ratios <- abs(est[sig]) / abs(true_effect)
  # Each carries its Monte Carlo error, NaN wherever the rate is, and Type M's needs two
  # significant replicates. The result stays a plain list, as code written for earlier versions
  # expects.
  type_s_ci <- if (usable) .type_s_interval(sum(wrong), n_sig, true_effect) else c(NaN, NaN)
  c(list(n_sims = n_sims, alpha = alpha, alternative = alternative),
    .rate_with_error(n_sig, n_sims, "power"),
    list(n_significant = n_sig, true_effect = true_effect, mean_estimate = mean(est),
         type_s = if (usable) mean(wrong) else NaN,
         type_s_lo = type_s_ci[1], type_s_hi = type_s_ci[2],
         type_m = if (usable) mean(ratios) else NaN,
         type_m_mcse = if (usable && n_sig >= 2L) .type_m_mcse(ratios) else NaN))
}

# The refusal for a design whose rows are correlated, byte-identical in the Python twin. It sends
# the user of either twin to R's power_mixed(), which fits every such design. The Python
# power_mixed() fits one only when it declares random effects.
.correlated_rows_message <- function(what) {
  sprintf(paste0(
    "power_design() in R and power() in Python t-test every row as an independent observation, ",
    "which is valid only when each subject contributes one row and no rows share a cluster. ",
    "This design has %s, so its rows are correlated and the t-test would overstate power and ",
    "understate Type M. Use power_mixed() in the R package, which fits the model the ",
    "specification implies."), what)
}

# Why the two-group backend cannot analyse a validated specification, or NULL when it can.
# A t-test of every row needs independent rows. It used to be applied to any Gaussian design with
# a two-level between factor, and under the null the test of 30 subjects crossed with 20 items was
# then significant in 105 of 200 replicates. An item unit, a within factor or an extra grouping
# factor is therefore refused, the first found being named. Predictors and a `subject` entry stay
# allowed, since with one row per subject they vary independently from row to row. A between
# factor on items needs an item unit, so the item test covers it. Mirrors _two_group_refusal()
# in power.py, message for message.
.two_group_refusal <- function(spec) {
  if (!identical(spec[["response"]][["family"]], "gaussian"))
    return("The power backend currently handles only the gaussian two-group design.")
  factors <- spec[["factors"]]
  between <- Filter(function(f) !is.null(f[["between"]]), factors)
  if (length(between) != 1L || length(between[[1L]][["levels"]]) != 2L)
    return("The power backend expects exactly one 2-level between factor.")
  # A factor can set `vary_within` alongside `between`, and it then varies within subjects.
  within <- Filter(function(f) !is.null(f[["vary_within"]]), factors)
  # The grouping factors simulate_design() draws: an `item` entry without an item unit groups
  # nothing, and with one the item unit is named first.
  extra <- setdiff(names(spec[["random"]]), c("subject", "item"))
  what <- if (!is.null(spec[["units"]][["item"]])) "an item unit"
          else if (length(within)) sprintf("the within factor '%s'", within[[1L]][["name"]])
          else if (length(extra)) sprintf("the grouping factor '%s'", extra[[1L]])
          else NULL
  if (is.null(what)) NULL else .correlated_rows_message(what)
}

# One Monte Carlo replicate of the two-group analysis. Kept at top level, out of a closure, so
# that only the arguments travel to PSOCK workers. Returns c(estimate, p). t.test() takes a
# one-sided p-value from the same t as the two-sided one, the upper tail for "greater" and the
# lower for "less", so it is half the two-sided p-value when t lies in the tested direction.
.power_design_rep <- function(i, spec, seeds, yname, fname, lev0, lev1,
                              alternative = "two.sided") {
  s <- spec; s[["seed"]] <- seeds[i]           # same seeds as the Python port
  d <- simulate_design(s, validate = FALSE)
  g0 <- d[[yname]][d[[fname]] == lev0]
  g1 <- d[[yname]][d[[fname]] == lev1]
  c(mean(g1) - mean(g0),
    stats::t.test(g1, g0, var.equal = TRUE, alternative = alternative)$p.value)
}
