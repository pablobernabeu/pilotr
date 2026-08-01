# Precision-based design analysis against a Region of Practical Equivalence (ROPE), together
# with an N-sweep to find the minimum analysable sample size. This is a fast frequentist
# analogue of a Bayesian HDI-versus-ROPE design analysis. Across Monte Carlo replicates we
# record, for each focal fixed effect, the 95% CI and whether it falls determinately outside
# the ROPE (a practically meaningful effect) or entirely inside it (practical equivalence).

#' Precision and ROPE design analysis at a fixed sample size
#'
#' A fast frequentist analogue of a Bayesian highest-density-interval-versus-ROPE design
#' analysis. Across Monte Carlo replicates, fit the model and record, for each focal fixed
#' effect, whether its 95% confidence interval falls entirely outside a region of practical
#' equivalence (a practically meaningful effect) or entirely inside it (practical equivalence
#' to zero), along with the expected interval width. Requires the `lme4` package.
#'
#' @details
#' The interval is a Wald approximation: the estimate plus or minus 1.96 standard errors
#' from the model's variance-covariance matrix. This fixed-z interval is chosen for speed and
#' for comparability across replicates; in small samples it is somewhat narrower than a
#' Satterthwaite t interval, so `p_meaningful` and `mean_ci_width` are slightly optimistic
#' at small sample sizes.
#'
#' @param spec A design specification (path or list).
#' @param focal The focal effects. `NULL`, the default, analyses every coefficient in the
#'   specification and takes the true values from it. A named numeric vector maps coefficient names
#'   to their true values, and a character vector names them without their true values. Interaction
#'   effects follow the model's column naming, so a specification key `a:b` is the focal name `a_b`.
#' @param formula Optional `lme4` formula; if `NULL` it is derived from the specification via
#'   [model_formula()].
#' @param prep Optional function mapping a simulated data set to the modelling data; if `NULL`
#'   it is derived via [model_data()], which log-transforms the outcome and
#'   builds the contrast and interaction columns, so focal names follow the
#'   auto-formula (interactions written as `a_b`).
#' @param rope Half-width of the region of practical equivalence; an effect with
#'   `abs(beta) < rope` is treated as practically equivalent to zero. Set it clearly narrower
#'   than the smallest effect worth detecting, because the probability of a determinate
#'   meaningful decision about an effect no larger than `rope` cannot rise above 0.5 however
#'   large the sample.
#' @param n_sims Number of Monte Carlo replicates. `p_meaningful` and `p_equivalent` are
#'   proportions over the converged replicates, so they carry a Monte Carlo standard error of
#'   about `sqrt(p * (1 - p) / n_sims)` and move in coarse steps when `n_sims` is small. At
#'   least 200 replicates are advisable for real planning.
#' @param workers Number of local worker processes over which to spread the replicates.
#'   The default of 1 runs serially. Because every replicate seeds the shared RNG from its
#'   own index, any worker count returns results identical to a serial run.
#' @return A data frame with one row per focal effect and columns `param`, `true`,
#'   `mean_ci_width`, `p_meaningful`, `p_equivalent`, `n_attempted`, `n_returned`,
#'   `n_converged`, `n_singular`, and `n_warning`. The interval behind
#'   `mean_ci_width` and the ROPE decisions is the Wald approximation described in
#'   Details.
#'
#'   The decision proportions are taken over `n_returned`, the replicates that
#'   produced an estimate. The remaining counts separate the fit outcomes, because a
#'   fit can return a usable estimate while still being boundary-singular or
#'   carrying a convergence warning: `n_converged` counts replicates with neither,
#'   `n_singular` those where `lme4::isSingular()` was true, and `n_warning` those
#'   with a warning or optimiser convergence message. Singular and warning fits are
#'   retained, since their fixed-effect estimates remain interpretable and
#'   discarding them would bias the result. A large `n_singular` means the model
#'   being fitted is richer than the design can support at that sample size, which
#'   is common in crossed designs (Bates et al., 2015; Matuschek et al., 2017).
#' @references Bates, D., Kliegl, R., Vasishth, S. and Baayen, H. (2015). Parsimonious mixed
#'   models. \emph{arXiv}. \doi{10.48550/arXiv.1506.04967}
#'
#'   Matuschek, H., Kliegl, R., Vasishth, S., Baayen, H. and Bates, D. (2017). Balancing
#'   Type I error and power in linear mixed models. \emph{Journal of Memory and Language},
#'   94, 305-315. \doi{10.1016/j.jml.2017.01.001}
#' @examples
#' \donttest{
#' if (requireNamespace("lme4", quietly = TRUE)) {
#'   spec <- build_spec(list(name = "pr", seed = 1, design_kind = "within",
#'     include_items = TRUE, n_subject = 12, n_item = 12, factor_name = "cond",
#'     lev1 = "a", lev2 = "b", intercept = 6, effect = 0.05,
#'     subj_int_sd = 0.12, subj_slope_sd = 0.04, subj_corr = 0.2,
#'     item_int_sd = 0.08, item_slope_sd = 0.02, item_corr = -0.1,
#'     family = "shifted_lognormal", resp_name = "", sigma = 0.3, shift = 200))
#'   # n_sims is small so the example runs quickly. Use 200 or more for real planning.
#'   precision_design(spec, focal = c(effect = 0.05), rope = 0.02, n_sims = 10)
#' }
#' }
#' @export
precision_design <- function(spec, focal = NULL, formula = NULL, prep = NULL, rope = 0.05,
                             n_sims = 100, workers = 1) {
  if (!requireNamespace("lme4", quietly = TRUE))
    stop("precision_design() requires the 'lme4' package; please install it.", call. = FALSE)
  spec <- .as_spec(spec)
  workers <- .check_workers(workers)
  cl <- NULL
  if (workers > 1L) {
    cl <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }
  .precision_design_impl(spec, focal, formula = formula, prep = prep, rope = rope,
                         n_sims = n_sims, cl = cl)
}

# The replicate loop behind precision_design(), taking an optional PSOCK cluster so that
# precision_curve() can create one cluster and reuse it across grid points.
.precision_design_impl <- function(spec, focal, formula, prep, rope, n_sims, cl = NULL) {
  if (is.null(formula)) formula <- model_formula(spec)
  if (is.null(prep)) prep <- .default_prep(spec)
  fo <- .resolve_focal(focal, spec)
  fnames <- fo$names
  if (!length(fnames))
    stop("this specification has no fixed coefficients to analyse; name the focal effects explicitly",
         call. = FALSE)

  seeds <- replicate_seeds(spec$seed, n_sims)
  # The shared replicate loop, with the cheaper fitter: a precision analysis needs estimates and
  # standard errors, and lmerTest's Satterthwaite p-values cost noticeably more than the plain fit.
  res <- .p_lapply(seq_len(n_sims), .design_rep, cl = cl, spec = spec, seeds = seeds,
                   prep = prep, formula = formula, fnames = fnames, test = FALSE)

  n_returned <- sum(vapply(res, function(r) isTRUE(r$fitted), logical(1)))
  coef_names <- NULL
  for (r in res) if (isTRUE(r$fitted)) { coef_names <- r$coef_names; break }

  seen <- setNames(integer(length(fnames)), fnames)
  blank <- setNames(rep(NA_real_, length(fnames)), fnames)
  width <- blank; p_out <- blank; out_mcse <- blank; out_lo <- blank; out_hi <- blank
  p_ins <- blank; ins_mcse <- blank; ins_lo <- blank; ins_hi <- blank

  for (f in fnames) {
    est <- vapply(res, function(r) r$est[[f]], numeric(1))
    se <- vapply(res, function(r) r$se[[f]], numeric(1))
    ok <- which(!is.na(est) & !is.na(se))
    seen[f] <- length(ok)
    if (!length(ok)) next
    lo <- est[ok] - .Z95 * se[ok]; hi <- est[ok] + .Z95 * se[ok]
    width[f] <- mean(hi - lo)
    n_out <- sum(lo > rope | hi < -rope)
    n_ins <- sum(lo > -rope & hi < rope)
    ro <- .rate_with_error(n_out, length(ok), "p")
    ri <- .rate_with_error(n_ins, length(ok), "p")
    p_out[f] <- ro$p; out_mcse[f] <- ro$p_mcse; out_lo[f] <- ro$p_lo; out_hi[f] <- ro$p_hi
    p_ins[f] <- ri$p; ins_mcse[f] <- ri$p_mcse; ins_lo[f] <- ri$p_lo; ins_hi[f] <- ri$p_hi
  }
  .warn_no_fits(res, n_returned, formula)
  .warn_absent_focal(seen, n_returned, coef_names)

  cnt <- .fit_counts(res, n_sims, n_returned)
  data.frame(
    param = fnames, true = unname(fo$true[fnames]),
    mean_ci_width = unname(width),
    p_meaningful = unname(p_out), p_meaningful_mcse = unname(out_mcse),
    p_meaningful_lo = unname(out_lo), p_meaningful_hi = unname(out_hi),
    p_equivalent = unname(p_ins), p_equivalent_mcse = unname(ins_mcse),
    p_equivalent_lo = unname(ins_lo), p_equivalent_hi = unname(ins_hi),
    n_attempted = cnt$n_attempted, n_returned = cnt$n_returned,
    n_converged = cnt$n_converged, n_singular = cnt$n_singular, n_warning = cnt$n_warning,
    row.names = NULL)
}

# The default data-preparation function, built in its own small frame so that only the
# spec (not the caller's whole frame, cluster included) travels to PSOCK workers.
.default_prep <- function(spec) function(d) model_data(spec, d)

#' Precision and ROPE curve over sample size
#'
#' Sweep the number of subjects and report the ROPE decision probabilities at each size, to
#' identify the minimum analysable \emph{N} at which a focal effect reaches a determinate decision
#' with a target probability (for example 0.90). Calls
#' [precision_design()] and so requires the `lme4` package.
#'
#' @details
#' A thin wrapper around [sweep_spec()] over `units$subject$n`, kept because a sample-size curve is
#' the sweep users want most often. For any other axis, call [sweep_spec()] directly; effect size is
#' the axis a design analysis most often needs next, and [design_conditions()] builds the
#' coefficient overrides for it.
#'
#' @param spec A design specification (path or list).
#' @param focal The focal effects, as in [precision_design()]. `NULL` uses every coefficient in the
#'   specification.
#' @param subject_ns A numeric vector of subject counts to evaluate.
#' @param formula Optional `lme4` formula; if `NULL` it is derived via
#'   [model_formula()].
#' @param prep Optional data-preparation function; if `NULL` it is derived via
#'   [model_data()].
#' @param rope Half-width of the region of practical equivalence. Set it clearly narrower than
#'   the smallest effect worth detecting, because the probability of a determinate meaningful
#'   decision about an effect no larger than `rope` cannot rise above 0.5 however large the
#'   sample, so the curve would fall with \emph{N} rather than rise.
#' @param n_sims Number of Monte Carlo replicates per sample size. `p_meaningful` and
#'   `p_equivalent` are proportions over the replicates that produced an estimate, so each is
#'   reported with its Monte Carlo standard error and Wilson interval. The default of 60 gives a
#'   standard error of 0.065 at a rate of 0.5, which is too coarse to support a claim about a
#'   design; raise it to at least 200 for planning.
#' @param workers Number of local worker processes over which to spread the replicates at
#'   each sample size. The default of 1 runs serially, and any worker count returns results
#'   identical to a serial run.
#' @return A data frame with one row per focal effect and sample size, adding an `n_subject`
#'   column to the columns returned by [precision_design()], including the Monte Carlo standard
#'   errors, the Wilson interval bounds, and the `n_returned`, `n_converged`, `n_singular` and
#'   `n_warning` fit counts.
#' @seealso [sweep_spec()] for any other axis, and [design_conditions()] for effect-size grids.
#' @examples
#' \donttest{
#' if (requireNamespace("lme4", quietly = TRUE)) {
#'   spec <- build_spec(list(name = "pr", seed = 1, design_kind = "within",
#'     include_items = TRUE, n_subject = 12, n_item = 12, factor_name = "cond",
#'     lev1 = "a", lev2 = "b", intercept = 6, effect = 0.05,
#'     subj_int_sd = 0.12, subj_slope_sd = 0.04, subj_corr = 0.2,
#'     item_int_sd = 0.08, item_slope_sd = 0.02, item_corr = -0.1,
#'     family = "shifted_lognormal", resp_name = "", sigma = 0.3, shift = 200))
#'   # n_sims is small so the example runs quickly. Use 200 or more for real planning.
#'   precision_curve(spec, focal = c(effect = 0.05), subject_ns = c(12, 18), rope = 0.02,
#'                   n_sims = 8)
#' }
#' }
#' @export
precision_curve <- function(spec, focal = NULL, subject_ns, formula = NULL, prep = NULL,
                            rope = 0.05, n_sims = 60, workers = 1) {
  if (!requireNamespace("lme4", quietly = TRUE))
    stop("precision_curve() requires the 'lme4' package; please install it.", call. = FALSE)
  sweep_spec(spec, "units$subject$n", subject_ns, precision_design,
             focal = focal, formula = formula, prep = prep, rope = rope,
             n_sims = n_sims, workers = workers, .name = "n_subject")
}
