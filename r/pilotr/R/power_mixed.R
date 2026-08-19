# Simulation-based power for mixed-effects designs, fit with lme4 / lmerTest.
# This capability distinguishes pilotr from the prototype (Bernabeu & Lynott, 2024,
# doi:10.5281/zenodo.10615953). It covers the territory of simr (Green & MacLeod, 2016)
# and mixedpower (Kumle et al., 2021), driven here by the portable design spec.
#
# Until 0.3 this backend accepted only a design with exactly one within-unit factor and an item
# unit, and fitted a formula written into the source rather than derived from the specification, so
# a design with two factors, with continuous predictors, or with a smaller random-effects structure
# than the maximal one was either refused outright or analysed under a model it had not described.
# It now runs the same replicate loop as precision_design(), over the model the specification
# implies.

#' Simulation-based power and design analysis for a mixed-effects design
#'
#' For each replicate, simulate from the ground-truth specification, fit the model the
#' specification implies with `lmerTest`, and test each focal fixed effect using Satterthwaite
#' p-values. Reports power together with the Type S and Type M errors of Gelman and Carlin (2014).
#' Requires the `lme4` and `lmerTest` packages.
#'
#' @details
#' `power_mixed()` is not a wrapper around an existing power package: it runs pilotr's own
#' simulation loop over the portable design specification. It covers territory pioneered by
#' `simr` (Green and MacLeod, 2016) and `mixedpower` (Kumle, Vo and Draschkow, 2021), to
#' which it is indebted. pilotr differs in being driven by the portable cross-language
#' specification, in reporting the Type S and Type M design-analysis errors alongside
#' power, and in parallelising its replicates through the `workers` argument.
#'
#' The analysis model comes from the specification rather than from this function. Before 0.3 the
#' formula was written into the source as a maximal crossed structure, so a design declaring
#' uncorrelated slopes, or no slopes at all, was nonetheless analysed as though it had them, and a
#' design with more than one factor was refused. The formula now comes from [model_formula()] and
#' the data from [model_data()], so the analysis matches the process that generated the data. Both
#' can still be given directly, which is what to do when a deliberately different analysis model is
#' the point, as when checking how a misspecified model behaves.
#'
#' Every reported rate carries its Monte Carlo standard error and a Wilson interval, because a
#' proportion over a finite number of replicates is an estimate rather than a fact. At the default
#' 100 replicates a power near 0.5 has a standard error of 0.05.
#'
#' @param spec A design specification (path or list).
#' @param focal The fixed effects to test. `NULL`, the default, tests every coefficient in the
#'   specification and takes the true values from it. A character vector names the effects and
#'   leaves the true values unknown, which suppresses Type S and Type M. A named numeric vector
#'   gives both, which is how to test against a value other than the one simulated. Interaction
#'   effects follow the model's column naming, so a specification key `a:b` is the focal name
#'   `a_b`.
#' @param formula Optional `lme4` formula; if `NULL` it is derived from the specification via
#'   [model_formula()].
#' @param prep Optional function mapping a simulated data set to the modelling data; if `NULL` it
#'   is derived via [model_data()].
#' @param n_sims Number of Monte Carlo replicates. A power estimate carries a Monte Carlo
#'   standard error of about `sqrt(p * (1 - p) / n_sims)`, and `type_s` and
#'   `type_m` average over the significant replicates alone, so they settle more
#'   slowly still. At least 200 replicates are advisable for study planning.
#' @param alpha Two-sided significance level.
#' @param workers Number of local worker processes over which to spread the replicates.
#'   The default of 1 runs serially. Because the replicate seeds are derived once from the
#'   specification's seed, any worker count returns results identical to a serial run. The
#'   mixed-model fits dominate the cost, so the speed-up is close to linear in the number of cores.
#' @return An object of class `pilotr_power`, a list whose per-run elements are `n_sims`,
#'   `alpha`, `n_attempted`, `n_returned`, `n_converged`, `n_singular` and `n_warning`, and whose
#'   per-effect elements are vectors named by focal effect: `power`, `power_mcse`, `power_lo`,
#'   `power_hi`, `n_significant`, `true_effect`, `mean_estimate`, `type_s` and `type_m`. With a
#'   single focal effect each of those has length one, so `result$power` reads as it always has.
#'
#'   `power` is the proportion of significant results among the replicates that returned an
#'   estimate for that effect, not among `n_sims`. The counts report the fit outcomes separately,
#'   because a fit can return a usable estimate while still being boundary-singular or carrying a
#'   convergence warning: `n_returned` counts replicates that yielded a fit, `n_converged` those
#'   that did so with neither a warning nor a singular fit, `n_singular` those where
#'   `lme4::isSingular()` was true, and `n_warning` those with a warning or optimiser convergence
#'   message. Singular and warning fits are retained in `power`, since their fixed-effect estimates
#'   remain interpretable and discarding them would bias the result: singularity is not independent
#'   of the variance estimates that produce it. A large `n_singular` means the model being fitted is
#'   richer than the design can support at that sample size, which is common in crossed designs
#'   (Bates et al., 2015; Matuschek et al., 2017), and is worth reporting alongside the power.
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
#'
#'   Bates, D., Kliegl, R., Vasishth, S. and Baayen, H. (2015). Parsimonious mixed models.
#'   \emph{arXiv}. \doi{10.48550/arXiv.1506.04967}
#'
#'   Matuschek, H., Kliegl, R., Vasishth, S., Baayen, H. and Bates, D. (2017). Balancing
#'   Type I error and power in linear mixed models. \emph{Journal of Memory and Language},
#'   94, 305-315. \doi{10.1016/j.jml.2017.01.001}
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
#'   # n_sims is small so the example runs quickly. Use 200 or more for real planning.
#'   power_mixed(spec, n_sims = 10)
#' }
#' }
#' @seealso [precision_design()] for the interval-width and ROPE analogue, and [sweep_spec()] to
#'   run this over a grid of sample sizes or effect sizes.
#' @export
power_mixed <- function(spec, focal = NULL, formula = NULL, prep = NULL,
                        n_sims = 100, alpha = 0.05, workers = 1) {
  if (!requireNamespace("lme4", quietly = TRUE) ||
      !requireNamespace("lmerTest", quietly = TRUE))
    stop("power_mixed() requires the 'lme4' and 'lmerTest' packages; please install them.",
         call. = FALSE)
  spec <- .as_spec(spec)
  workers <- .check_workers(workers)
  cl <- NULL
  if (workers > 1L) {
    cl <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }
  .power_mixed_impl(spec, focal = focal, formula = formula, prep = prep,
                    n_sims = n_sims, alpha = alpha, cl = cl)
}

# The replicate loop behind power_mixed(), taking an optional PSOCK cluster so that sweep
# functions can create one cluster and reuse it across grid points.
.power_mixed_impl <- function(spec, focal, formula, prep, n_sims, alpha, cl = NULL) {
  if (is.null(formula)) formula <- model_formula(spec)
  if (is.null(prep)) prep <- .default_prep(spec)
  fo <- .resolve_focal(focal, spec)
  fnames <- fo$names
  if (!length(fnames))
    stop("this specification has no fixed coefficients to test; name the focal effects explicitly",
         call. = FALSE)

  seeds <- replicate_seeds(spec$seed, n_sims)
  res <- .p_lapply(seq_len(n_sims), .design_rep, cl = cl, spec = spec, seeds = seeds,
                   prep = prep, formula = formula, fnames = fnames, test = TRUE)

  n_returned <- sum(vapply(res, function(r) isTRUE(r$fitted), logical(1)))
  coef_names <- NULL
  for (r in res) if (isTRUE(r$fitted)) { coef_names <- r$coef_names; break }

  seen <- stats::setNames(integer(length(fnames)), fnames)
  blank <- stats::setNames(rep(NA_real_, length(fnames)), fnames)
  power <- mcse <- lo <- hi <- type_s <- type_m <- mean_est <- blank
  n_sig <- stats::setNames(integer(length(fnames)), fnames)

  for (f in fnames) {
    est <- vapply(res, function(r) r$est[[f]], numeric(1))
    pv <- vapply(res, function(r) r$p[[f]], numeric(1))
    ok <- which(!is.na(pv))
    seen[f] <- length(ok)
    if (!length(ok)) next
    sig <- ok[pv[ok] < alpha]
    rate <- .rate_with_error(length(sig), length(ok), "power")
    power[f] <- rate$power; mcse[f] <- rate$power_mcse
    lo[f] <- rate$power_lo; hi[f] <- rate$power_hi
    n_sig[f] <- length(sig)
    mean_est[f] <- mean(est[ok])
    beta <- fo$true[[f]]
    # Type S and Type M are defined relative to a true value, and Type M divides by it, so both
    # stay NA when the true effect is unknown or zero rather than reporting an infinity.
    if (length(sig) && !is.na(beta) && beta != 0) {
      type_s[f] <- mean((est[sig] > 0) != (beta > 0))
      type_m[f] <- mean(abs(est[sig]) / abs(beta))
    }
  }
  .warn_no_fits(res, n_returned, formula)
  .warn_absent_focal(seen, n_returned, coef_names)

  out <- c(list(n_sims = n_sims, alpha = alpha),
           .fit_counts(res, n_sims, n_returned),
           list(power = power, power_mcse = mcse, power_lo = lo, power_hi = hi,
                n_significant = n_sig, true_effect = fo$true[fnames],
                mean_estimate = mean_est, type_s = type_s, type_m = type_m))
  class(out) <- c("pilotr_power", "list")
  out
}

#' Print a simulation-based power result
#'
#' Shows the fit accounting and then one row per focal effect, with each power estimate beside its
#' Monte Carlo standard error and Wilson interval, so that the precision of the estimate is as
#' visible as the estimate.
#'
#' @param x A `pilotr_power` object, as returned by [power_mixed()].
#' @param digits Number of significant digits for the reported rates.
#' @param ... Ignored, present for consistency with the generic.
#' @return `x`, invisibly.
#' @export
print.pilotr_power <- function(x, digits = 3, ...) {
  cat(sprintf("Simulation-based power over %d replicates (alpha = %g)\n", x$n_sims, x$alpha))
  cat(sprintf("  fits: %d attempted, %d returned, %d converged cleanly, %d singular, %d with warnings\n",
              x$n_attempted, x$n_returned, x$n_converged, x$n_singular, x$n_warning))
  if (x$n_singular > 0.25 * max(x$n_attempted, 1))
    cat("  note: many fits were boundary-singular, so this model is richer than the design supports\n")
  tab <- data.frame(
    effect = names(x$power),
    true = unname(x$true_effect),
    power = round(unname(x$power), digits),
    mcse = round(unname(x$power_mcse), digits),
    ci95 = ifelse(is.na(x$power_lo), NA_character_,
                  sprintf("[%.3f, %.3f]", unname(x$power_lo), unname(x$power_hi))),
    n_sig = unname(x$n_significant),
    type_s = round(unname(x$type_s), digits),
    type_m = round(unname(x$type_m), digits),
    row.names = NULL, stringsAsFactors = FALSE)
  print(tab, row.names = FALSE)
  invisible(x)
}

#' Power curve over sample size for a mixed-effects design
#'
#' Sweep the number of subjects and compute mixed-effects power at each. Pass the result to
#' [target_n()] for the sample size at which power crosses a target, with an interval on it. A
#' thin wrapper around [sweep_spec()] over `units$subject$n`, kept because a sample-size curve is
#' the sweep users want most often. Requires the `lme4` and `lmerTest` packages.
#'
#' @details
#' Like [power_mixed()], this runs pilotr's own simulation loop over the portable design
#' specification rather than wrapping an existing package, and differs from `simr` (Green and
#' MacLeod, 2016) and `mixedpower` (Kumle, Vo and Draschkow, 2021) in being driven by that
#' specification, in reporting Type S and Type M errors, and in built-in parallelisation: with
#' `workers > 1` a single worker pool is created once and reused across all sample sizes.
#'
#' For any axis other than sample size, call [sweep_spec()] directly. Effect size is the axis a
#' design analysis most often needs after sample size, and [design_conditions()] builds the
#' coefficient overrides for it.
#'
#' The curve is the input to [target_n()], which is where the sample size a preregistration
#' quotes should come from. Reading the crossing off a plot instead judges points whose Monte
#' Carlo intervals overlap, and reports the answer without the interval that goes with it.
#'
#' @param spec A design specification (path or list).
#' @param subject_ns A numeric vector of subject counts to evaluate.
#' @param focal The fixed effects to test, as in [power_mixed()].
#' @param n_sims Number of Monte Carlo replicates per point. A power estimate carries a Monte
#'   Carlo standard error of about `sqrt(p * (1 - p) / n_sims)`, reported alongside it as
#'   `power_mcse`. The default of 60 gives a standard error of 0.065 at a power of 0.5, which is
#'   too coarse to support a claim about a design; raise it to at least 200 for planning.
#' @param alpha Two-sided significance level.
#' @param workers Number of local worker processes over which to spread the replicates at
#'   each grid point. The default of 1 runs serially, and any worker count returns results
#'   identical to a serial run.
#' @return A data frame with one row per sample size and focal effect, with columns `n_subject`,
#'   `effect`, `true`, `power`, `power_mcse`, `power_lo`, `power_hi`, `n_significant`, `type_s`,
#'   `type_m`, and the `n_attempted`, `n_returned`, `n_converged`, `n_singular` and `n_warning` fit
#'   counts. `n_singular` typically falls as the sample size rises, so reading it down the sweep
#'   shows where the model becomes supportable.
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
#'   # n_sims is small so the example runs quickly. Use 200 or more for real planning.
#'   power_curve_mixed(spec, subject_ns = c(12, 18), n_sims = 8)
#' }
#' }
#' @seealso [target_n()] to solve the returned curve for a sample size, [sweep_spec()] for any
#'   other axis, and [design_conditions()] for effect-size grids.
#' @export
power_curve_mixed <- function(spec, subject_ns, focal = NULL, n_sims = 60, alpha = 0.05,
                              workers = 1) {
  if (!requireNamespace("lme4", quietly = TRUE) ||
      !requireNamespace("lmerTest", quietly = TRUE))
    stop("power_curve_mixed() requires the 'lme4' and 'lmerTest' packages; please install them.",
         call. = FALSE)
  sweep_spec(spec, "units$subject$n", subject_ns, power_mixed,
             focal = focal, n_sims = n_sims, alpha = alpha, workers = workers,
             .name = "n_subject")
}

# Flatten a pilotr_power object into one data-frame row per focal effect, for the sweeps.
.as_power_frame <- function(x) {
  data.frame(
    effect = names(x$power), true = unname(x$true_effect),
    power = unname(x$power), power_mcse = unname(x$power_mcse),
    power_lo = unname(x$power_lo), power_hi = unname(x$power_hi),
    n_significant = unname(x$n_significant),
    type_s = unname(x$type_s), type_m = unname(x$type_m),
    n_attempted = x$n_attempted, n_returned = x$n_returned, n_converged = x$n_converged,
    n_singular = x$n_singular, n_warning = x$n_warning,
    row.names = NULL, stringsAsFactors = FALSE)
}
