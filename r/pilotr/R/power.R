# Simulation-based power and design analysis (Type S and Type M), mirroring power.py.

#' Simulation-based power and design analysis for a two-group Gaussian design
#'
#' Estimate power by repeatedly simulating from the specification and applying a two-sample
#' t-test, alongside the Type S (sign) and Type M (magnitude) design-analysis errors of
#' Gelman and Carlin (2014).
#'
#' @param spec A design specification (path or list) for a two-group Gaussian design.
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
  if (spec$response$family != "gaussian")
    stop("The power backend currently handles only the gaussian two-group design.")
  between <- Filter(function(f) !is.null(f$between), spec$factors)
  if (length(between) != 1 || length(between[[1]]$levels) != 2)
    stop("The power backend expects exactly one 2-level between factor.")
  f <- between[[1]]; fname <- f$name; lev0 <- f$levels[1]; lev1 <- f$levels[2]
  col <- names(f$contrasts)[1]; vals <- f$contrasts[[col]]
  true_effect <- spec$fixed$coefficients[[col]] * (vals[2] - vals[1])
  yname <- spec$response$name
  seeds <- replicate_seeds(spec$seed, n_sims)

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

# One Monte Carlo replicate of the two-group analysis. Kept at top level, out of a closure, so
# that only the arguments travel to PSOCK workers. Returns c(estimate, p).
.power_design_rep <- function(i, spec, seeds, yname, fname, lev0, lev1) {
  s <- spec; s$seed <- seeds[i]                # same seeds as the Python port
  d <- simulate_design(s, validate = FALSE)
  g0 <- d[[yname]][d[[fname]] == lev0]
  g1 <- d[[yname]][d[[fname]] == lev1]
  c(mean(g1) - mean(g0), stats::t.test(g1, g0, var.equal = TRUE)$p.value)
}
