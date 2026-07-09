# Simulation-based power for crossed mixed-effects designs, fit with lme4 / lmerTest.
# This capability distinguishes pilotr from the prototype (Bernabeu & Lynott, 2024,
# doi:10.5281/zenodo.10615953). It covers the territory of simr (Green & MacLeod, 2016)
# and mixedpower (Kumle et al., 2021), driven here by the portable design spec. v0.1
# handles one within-subject/within-item factor (single contrast column) with a gaussian
# or shifted_lognormal response.

#' Simulation-based power and design analysis for a crossed mixed-effects design
#'
#' For each replicate, simulate from the ground-truth specification, fit the maximal model
#' `y ~ cond + (1 + cond | subject) + (1 + cond | item)` with `lmerTest`, and test the fixed
#' effect using Satterthwaite p-values. Reports power together with the Type S and Type M
#' errors of Gelman and Carlin (2014). Requires the `lme4` and `lmerTest` packages.
#'
#' @details
#' `power_mixed()` is not a wrapper around an existing power package: it runs pilotr's own
#' simulation loop over the portable design specification. It covers territory pioneered by
#' `simr` (Green and MacLeod, 2016) and `mixedpower` (Kumle, Vo and Draschkow, 2021), to
#' which it is indebted. pilotr differs in being driven by the portable cross-language
#' specification, so that R and Python produce bit-identical data, in reporting the Type S
#' and Type M design-analysis errors alongside power, and in parallelising its replicates
#' through the `workers` argument.
#'
#' @param spec A design specification (path or list) with exactly one within-unit factor.
#' @param n_sims Number of Monte Carlo replicates.
#' @param alpha Two-sided significance level.
#' @param workers Number of local worker processes over which to spread the replicates.
#'   The default of 1 runs serially. Because every replicate seeds the shared RNG from its
#'   own index, any worker count returns results identical to a serial run. The mixed-model
#'   fits dominate the cost, so the speed-up is close to linear in the number of cores.
#' @return A list with elements `n_sims`, `n_converged`, `alpha`, `power`, `n_significant`,
#'   `true_effect`, `mean_estimate`, `type_s`, and `type_m`.
#' @references Gelman, A. and Carlin, J. (2014). Beyond power calculations: Assessing Type S
#'   (sign) and Type M (magnitude) errors. \emph{Perspectives on Psychological Science},
#'   9(6), 641-651. \doi{10.1177/1745691614551642}
#'
#'   Green, P. and MacLeod, C. J. (2016). SIMR: An R package for power analysis of
#'   generalized linear mixed models by simulation. \emph{Methods in Ecology and Evolution},
#'   7(4), 493-498. \doi{10.1111/2041-210x.12504}
#'
#'   Kumle, L., Vo, M. L.-H. and Draschkow, D. (2021). Estimating power in (generalized)
#'   linear mixed models: An open introduction and tutorial in R. \emph{Behavior Research
#'   Methods}, 53, 2528-2543. \doi{10.3758/s13428-021-01546-0}
#' @examples
#' \donttest{
#' if (requireNamespace("lme4", quietly = TRUE) &&
#'     requireNamespace("lmerTest", quietly = TRUE)) {
#'   spec <- build_spec(list(name = "p", seed = 1, design_kind = "within",
#'     include_items = TRUE, n_subject = 12, n_item = 12, factor_name = "cond",
#'     lev1 = "a", lev2 = "b", intercept = 6, effect = 0.05,
#'     subj_int_sd = 0.12, subj_slope_sd = 0.04, subj_corr = 0.2,
#'     item_int_sd = 0.08, item_slope_sd = 0.02, item_corr = -0.1,
#'     family = "shifted_lognormal", resp_name = "", sigma = 0.3, shift = 200))
#'   power_mixed(spec, n_sims = 10)
#' }
#' }
#' @export
power_mixed <- function(spec, n_sims = 100, alpha = 0.05, workers = 1) {
  if (!requireNamespace("lme4", quietly = TRUE) ||
      !requireNamespace("lmerTest", quietly = TRUE))
    stop("power_mixed() requires the 'lme4' and 'lmerTest' packages; please install them.",
         call. = FALSE)
  workers <- .check_workers(workers)
  cl <- NULL
  if (workers > 1L) {
    cl <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }
  .power_mixed_impl(spec, n_sims = n_sims, alpha = alpha, cl = cl)
}

# The replicate loop behind power_mixed(), taking an optional PSOCK cluster so that sweep
# functions can create one cluster and reuse it across grid points.
.power_mixed_impl <- function(spec, n_sims, alpha, cl = NULL) {
  if (is.character(spec)) spec <- load_spec(spec)
  within <- Filter(function(f) !is.null(f$vary_within), spec$factors)
  if (length(within) != 1) stop("v0.1 power_mixed expects exactly one within factor.")
  f <- within[[1]]; fname <- f$name
  col <- names(f$contrasts)[1]; cvals <- f$contrasts[[col]]; levs <- f$levels
  beta <- spec$fixed$coefficients[[col]]
  yname <- spec$response$name
  fam <- spec$response$family
  shift <- if (is.null(spec$response$shift)) 0 else spec$response$shift
  base <- spec$seed

  res <- .p_lapply(seq_len(n_sims), .power_mixed_rep, cl = cl, spec = spec, base = base,
                   fname = fname, levs = levs, cvals = cvals, yname = yname, fam = fam,
                   shift = shift)
  est <- vapply(res, `[[`, numeric(1), 1L)
  pv  <- vapply(res, `[[`, numeric(1), 2L)
  ok <- which(!is.na(pv))
  sig <- ok[pv[ok] < alpha]
  list(
    n_sims = n_sims, n_converged = length(ok), alpha = alpha,
    power = length(sig) / length(ok), n_significant = length(sig),
    true_effect = beta, mean_estimate = mean(est[ok]),
    type_s = if (length(sig)) mean((est[sig] > 0) != (beta > 0)) else NaN,
    type_m = if (length(sig)) mean(abs(est[sig]) / abs(beta)) else NaN
  )
}

# One mixed-model replicate. Kept at top level so that only the arguments travel to PSOCK
# workers. Returns c(estimate, p), or c(NA, NA) when the fit fails.
.power_mixed_rep <- function(i, spec, base, fname, levs, cvals, yname, fam, shift) {
  s <- spec; s$seed <- base + (i - 1)
  d <- simulate_design(s)
  d$.cond <- cvals[match(d[[fname]], levs)]
  d$.y <- if (fam == "shifted_lognormal") log(d[[yname]] - shift) else d[[yname]]
  fit <- tryCatch(
    suppressWarnings(suppressMessages(
      lmerTest::lmer(.y ~ .cond + (1 + .cond | subject) + (1 + .cond | item),
                     data = d, control = lme4::lmerControl(calc.derivs = FALSE)))),
    error = function(e) NULL)
  if (is.null(fit)) return(c(NA_real_, NA_real_))
  co <- summary(fit)$coefficients
  c(co[".cond", "Estimate"], co[".cond", "Pr(>|t|)"])
}

#' Power curve over sample size for a crossed mixed-effects design
#'
#' Sweep the number of subjects and compute mixed-effects power at each, to locate where
#' power crosses a target. Runs the same replicate loop as [power_mixed()] and so requires
#' the `lme4` and `lmerTest` packages.
#'
#' @details
#' Like [power_mixed()], this function runs pilotr's own simulation loop over the portable
#' design specification rather than wrapping an existing package. It covers territory
#' pioneered by `simr` (Green and MacLeod, 2016) and `mixedpower` (Kumle, Vo and Draschkow,
#' 2021), and differs in being driven by the portable cross-language specification with
#' bit-identical R and Python output, in reporting Type S and Type M errors, and in
#' built-in parallelisation: with `workers > 1` a single worker pool is created once and
#' reused across all sample sizes in the sweep.
#'
#' @param spec A design specification (path or list) with one within-unit factor.
#' @param subject_ns A numeric vector of subject counts to evaluate.
#' @param n_sims Number of Monte Carlo replicates per point.
#' @param alpha Two-sided significance level.
#' @param workers Number of local worker processes over which to spread the replicates at
#'   each grid point. The default of 1 runs serially, and any worker count returns results
#'   identical to a serial run.
#' @return A data frame with one row per sample size and columns `n_subject`, `power`, and
#'   `type_m`.
#' @references Green, P. and MacLeod, C. J. (2016). SIMR: An R package for power analysis
#'   of generalized linear mixed models by simulation. \emph{Methods in Ecology and
#'   Evolution}, 7(4), 493-498. \doi{10.1111/2041-210x.12504}
#'
#'   Kumle, L., Vo, M. L.-H. and Draschkow, D. (2021). Estimating power in (generalized)
#'   linear mixed models: An open introduction and tutorial in R. \emph{Behavior Research
#'   Methods}, 53, 2528-2543. \doi{10.3758/s13428-021-01546-0}
#' @examples
#' \donttest{
#' if (requireNamespace("lme4", quietly = TRUE) &&
#'     requireNamespace("lmerTest", quietly = TRUE)) {
#'   spec <- build_spec(list(name = "p", seed = 1, design_kind = "within",
#'     include_items = TRUE, n_subject = 12, n_item = 12, factor_name = "cond",
#'     lev1 = "a", lev2 = "b", intercept = 6, effect = 0.05,
#'     subj_int_sd = 0.12, subj_slope_sd = 0.04, subj_corr = 0.2,
#'     item_int_sd = 0.08, item_slope_sd = 0.02, item_corr = -0.1,
#'     family = "shifted_lognormal", resp_name = "", sigma = 0.3, shift = 200))
#'   power_curve_mixed(spec, subject_ns = c(12, 18), n_sims = 8)
#' }
#' }
#' @export
power_curve_mixed <- function(spec, subject_ns, n_sims = 60, alpha = 0.05, workers = 1) {
  if (!requireNamespace("lme4", quietly = TRUE) ||
      !requireNamespace("lmerTest", quietly = TRUE))
    stop("power_curve_mixed() requires the 'lme4' and 'lmerTest' packages; please install them.",
         call. = FALSE)
  if (is.character(spec)) spec <- load_spec(spec)
  workers <- .check_workers(workers)
  cl <- NULL
  if (workers > 1L) {
    cl <- parallel::makeCluster(workers)   # one cluster, reused across the whole sweep
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }
  parts <- lapply(subject_ns, function(n) {
    s <- spec; s$units$subject$n <- n
    r <- .power_mixed_impl(s, n_sims = n_sims, alpha = alpha, cl = cl)
    data.frame(n_subject = n, power = r$power, type_m = r$type_m)
  })
  do.call(rbind, parts)
}
