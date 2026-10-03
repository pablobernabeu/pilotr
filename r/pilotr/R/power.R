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
#' @param spec A design specification (path or list) for a two-group Gaussian design with one
#'   row per subject (see Details).
#' @param n_sims Number of Monte Carlo replicates. A power estimate carries a Monte Carlo
#'   standard error of about `sqrt(p * (1 - p) / n_sims)`, and `type_s` and
#'   `type_m` average over the significant replicates alone, so they settle more
#'   slowly still. At least 200 replicates are advisable for study planning.
#' @param alpha Two-sided significance level.
#' @param workers Number of local worker processes over which to spread the replicates.
#'   The default of 1 runs serially. Because every replicate takes its own seed from
#'   [replicate_seeds()], any worker count returns results identical to a serial run.
#' @return A list with elements `n_sims`, `alpha`, `power`, `n_significant`,
#'   `true_effect`, `mean_estimate`, `type_s` (sign-error rate among significant
#'   replicates), and `type_m` (mean exaggeration ratio among significant
#'   replicates). Both design-analysis quantities are `NaN` when no replicate reached
#'   significance and when the true effect is zero, as in the null condition
#'   [design_conditions()] produces: neither is defined without a true value to
#'   compare against, and Type M divides by it.
#' @references Gelman, A. and Carlin, J. (2014). Beyond power calculations: Assessing Type S
#'   (sign) and Type M (magnitude) errors. \emph{Perspectives on Psychological Science},
#'   9(6), 641-651. \doi{10.1177/1745691614551642}
#' @examples
#' spec <- build_spec(list(name = "d", seed = 1, design_kind = "between",
#'   factor_name = "group", lev1 = "a", lev2 = "b", n_subject = 64,
#'   intercept = 100, effect = 5, family = "gaussian", resp_name = "", sigma = 10))
#' # n_sims is small so the example runs quickly. Use 200 or more for real planning.
#' power_design(spec, n_sims = 50)
#' @export
power_design <- function(spec, n_sims = 1000, alpha = 0.05, workers = 1) {
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
  cl <- NULL
  if (workers > 1L) {
    cl <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }
  res <- .p_lapply(seq_len(n_sims), .power_design_rep, cl = cl, spec = spec, seeds = seeds,
                   yname = yname, fname = fname, lev0 = lev0, lev1 = lev1)
  est <- vapply(res, `[[`, numeric(1), 1L)
  pv  <- vapply(res, `[[`, numeric(1), 2L)
  sig <- which(pv < alpha)
  # Type S and Type M are defined relative to a true value, and Type M divides by it, so both
  # stay NaN when the true effect is zero. The alternative was an infinity, alongside a
  # sign-error rate that had quietly become "the estimate is positive". Same rule as
  # power_mixed(), and the null condition design_conditions() recommends is exactly this case.
  usable <- length(sig) > 0L && !is.na(true_effect) && true_effect != 0
  list(
    n_sims = n_sims, alpha = alpha,
    power = length(sig) / n_sims, n_significant = length(sig),
    true_effect = true_effect, mean_estimate = mean(est),
    type_s = if (usable) mean((est[sig] > 0) != (true_effect > 0)) else NaN,
    type_m = if (usable) mean(abs(est[sig]) / abs(true_effect)) else NaN
  )
}

# The refusal for a design whose rows are correlated, byte-identical in the Python twin. It sends
# the user of either twin to R's power_mixed(), since the Python power_mixed() fits a single
# within factor crossed with items and refuses every between-subjects design.
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
# that only the arguments travel to PSOCK workers. Returns c(estimate, p).
.power_design_rep <- function(i, spec, seeds, yname, fname, lev0, lev1) {
  s <- spec; s[["seed"]] <- seeds[i]           # same seeds as the Python port
  d <- simulate_design(s, validate = FALSE)
  g0 <- d[[yname]][d[[fname]] == lev0]
  g1 <- d[[yname]][d[[fname]] == lev1]
  c(mean(g1) - mean(g0), stats::t.test(g1, g0, var.equal = TRUE)$p.value)
}
