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
#' @param spec A design specification (path or list) with exactly one within-unit factor.
#' @param n_sims Number of Monte Carlo replicates.
#' @param alpha Two-sided significance level.
#' @return A list with elements `n_sims`, `n_converged`, `alpha`, `power`, `n_significant`,
#'   `true_effect`, `mean_estimate`, `type_s`, and `type_m`.
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
power_mixed <- function(spec, n_sims = 100, alpha = 0.05) {
  if (!requireNamespace("lme4", quietly = TRUE) ||
      !requireNamespace("lmerTest", quietly = TRUE))
    stop("power_mixed() requires the 'lme4' and 'lmerTest' packages; please install them.",
         call. = FALSE)
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

  est <- rep(NA_real_, n_sims); pv <- rep(NA_real_, n_sims)
  for (i in seq_len(n_sims)) {
    s <- spec; s$seed <- base + (i - 1)
    d <- simulate_design(s)
    d$.cond <- cvals[match(d[[fname]], levs)]
    d$.y <- if (fam == "shifted_lognormal") log(d[[yname]] - shift) else d[[yname]]
    fit <- tryCatch(
      suppressWarnings(suppressMessages(
        lmerTest::lmer(.y ~ .cond + (1 + .cond | subject) + (1 + .cond | item),
                       data = d, control = lme4::lmerControl(calc.derivs = FALSE)))),
      error = function(e) NULL)
    if (is.null(fit)) next
    co <- summary(fit)$coefficients
    est[i] <- co[".cond", "Estimate"]
    pv[i]  <- co[".cond", "Pr(>|t|)"]
  }
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

#' Power curve over sample size for a crossed mixed-effects design
#'
#' Sweep the number of subjects and compute mixed-effects power at each, to locate where
#' power crosses a target. Calls [power_mixed()] and so requires the `lme4` and `lmerTest`
#' packages.
#'
#' @param spec A design specification (path or list) with one within-unit factor.
#' @param subject_ns A numeric vector of subject counts to evaluate.
#' @param n_sims Number of Monte Carlo replicates per point.
#' @param alpha Two-sided significance level.
#' @return A data frame with one row per sample size and columns `n_subject`, `power`, and
#'   `type_m`.
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
power_curve_mixed <- function(spec, subject_ns, n_sims = 60, alpha = 0.05) {
  if (is.character(spec)) spec <- load_spec(spec)
  parts <- lapply(subject_ns, function(n) {
    s <- spec; s$units$subject$n <- n
    r <- power_mixed(s, n_sims = n_sims, alpha = alpha)
    data.frame(n_subject = n, power = r$power, type_m = r$type_m)
  })
  do.call(rbind, parts)
}
