# Precision-based design analysis against a Region of Practical Equivalence (ROPE), together
# with an N-sweep over sample size. This is a fast frequentist analogue of a Bayesian
# HDI-versus-ROPE design analysis. Across Monte Carlo replicates we record, for each focal
# fixed effect, the 95% CI and whether it falls determinately outside the ROPE (a practically
# meaningful effect) or entirely inside it (practical equivalence).
#
# The sweep reports the decision probability at each sample size. Turning that into the minimum
# analysable N is solve_curve()'s job, not the reader's: at the replicate counts these functions
# are usually run at, neighbouring points are not significantly different, so a crossing judged
# by eye carries an uncertainty nobody has written down.

#' Precision and ROPE design analysis at a fixed sample size
#'
#' A fast frequentist analogue of a Bayesian highest-density-interval-versus-ROPE design
#' analysis. Across Monte Carlo replicates, fit the model and record, for each focal fixed
#' effect, whether its 95% confidence interval falls entirely outside a region of practical
#' equivalence (a practically meaningful effect) or entirely inside it (practical equivalence),
#' along with the expected interval width. Requires the `lme4` package.
#'
#' @details
#' The interval is a Wald approximation: the estimate plus or minus 1.96 standard errors
#' from the model's variance-covariance matrix. This fixed-z interval is chosen for speed and
#' for comparability across replicates; in small samples it is somewhat narrower than a
#' Satterthwaite t interval, so `p_meaningful` and `mean_ci_width` are slightly optimistic
#' at small sample sizes.
#'
#' The fitter follows the response family, as in [power_mixed()]. Since no p-value is needed, a
#' `gaussian`, `lognormal`, `shifted_lognormal` or `exgaussian` response is fitted by
#' `lme4::lmer()`. A `bernoulli` or `poisson` response is fitted by `lme4::glmer()` with the
#' binomial or Poisson family, and a model with no random terms by `stats::lm()` or
#' `stats::glm()`. The interval is therefore on the scale of the specification's coefficients,
#' and so is the region of practical equivalence that `rope` sets. That is the log scale for the
#' lognormal families and `poisson`, and the logit scale for `bernoulli`. A `glmer()` fit has the
#' Wald interval that Bolker et al. (2009) recommend for a GLMM without overdispersion. With few
#' subjects or items, it is likewise narrower than it should be, since the test it inverts rejects
#' too often when the clusters are few (Li and Redden, 2015). The `fitter` column names the
#' fitter.
#'
#' pilotr has no frequentist model for an `ordinal` or `beta` response yet. Such a response is
#' fitted by the linear model on its own scale, where a region of practical equivalence set on
#' the logit scale of the coefficients means nothing. The function warns, and withholds the
#' decision probabilities with their Monte Carlo standard errors and bounds (`NA`).
#' `mean_ci_width` stays, as the width of an interval on the response scale.
#' [generate_design_analysis()] applies an interval-and-ROPE rule on the link scale, in a
#' Bayesian model.
#'
#' A `formula` given here is fitted as written, by `lme4::lmer()`, or by `stats::lm()` when it
#' has no random terms, whatever the family. With such a formula, or a `prep` that changes the
#' response, and `focal = NULL`, the true values still come from the specification and are on its
#' scale, which need not be the model's. A named numeric `focal` gives the true values on the
#' model's scale.
#'
#' `formula` given as a ladder of models simulates a preregistered simplification rule, as in
#' [power_mixed()]. Each replicate fits the rungs in turn and keeps the first fit that `fallback`
#' accepts, or the last rung's fit when it accepts none. The intervals and decisions are then those
#' of the rule, and one column per rung gives the share of the replicates analysed at it.
#'
#' @param spec A design specification (path or list).
#' @param focal The focal effects. `NULL`, the default, analyses every coefficient in the
#'   specification and takes the true values from it. A named numeric vector maps coefficient names
#'   to their true values, and a character vector names them without their true values. Interaction
#'   effects follow the model's column naming, so a specification key `a:b` is the focal name `a_b`.
#'   A list of coefficient sets, which [power_mixed()] tests jointly, is refused, since a set has no
#'   single estimate to build an interval around.
#' @param formula Optional `lme4` formula; if `NULL` it is derived from the specification via
#'   [model_formula()], and the response family chooses the fitter. A formula given here is
#'   fitted as a linear model whatever the family (see Details). A ladder of models simulates a
#'   simplification rule: a character vector of [model_formula()] structures, such as
#'   `c("maximal", "zcp", "intercepts")`, fitted as the derived model is, or a list of formulas,
#'   fitted as written, as in [power_mixed()].
#' @param fallback When a replicate moves on to the next rung of a ladder: `"singular"`, the
#'   default, after a fit that is boundary-singular or failed, and `"warning"` also after a fit with
#'   a warning or an optimiser convergence message. It has no effect on a single formula.
#' @param prep Optional function mapping a simulated data set to the modelling data; if `NULL`
#'   it is derived via [model_data()], which adds the analysis response `.y`, on the log scale
#'   for the lognormal families, and builds the contrast and interaction columns, so focal names
#'   follow the auto-formula (interactions written as `a_b`).
#' @param rope The region of practical equivalence, as a half-width or as an interval. A single
#'   non-negative number `r` is the region from `-r` to `r`, so that an effect with
#'   `abs(beta) < r` is treated as practically equivalent to zero. A pair `c(lower, upper)` with
#'   `lower < upper` is the region between those bounds, which need not be centred on zero, as
#'   for a prediction that an effect lies within a stated interval. Any other value is refused.
#'   Keep the effects worth detecting clearly outside the region, because the probability of a
#'   determinate meaningful decision about an effect inside it cannot rise above 0.5 however
#'   large the sample.
#' @param n_sims Number of Monte Carlo replicates. `p_meaningful` and `p_equivalent` are
#'   proportions over the replicates that returned an estimate (`n_returned`), so they carry
#'   a Monte Carlo standard error of about `sqrt(p * (1 - p) / n_sims)` and move in coarse
#'   steps when `n_sims` is small. At least 200 replicates are advisable for real planning.
#' @param workers Number of local worker processes over which to spread the replicates.
#'   The default of 1 runs serially. Because every replicate takes its own seed from
#'   [replicate_seeds()], any worker count returns results identical to a serial run.
#' @param cl A cluster made by [parallel::makeCluster()] to run the replicates on, in place of one
#'   started for this call, in which case `workers` is ignored. The cluster is left running for its
#'   owner to stop, and through [sweep_spec()] one cluster serves every grid point.
#' @return A data frame with one row per focal effect and columns `param`, `true`,
#'   `mean_ci_width`, `p_meaningful`, `p_meaningful_mcse`, `p_meaningful_lo`,
#'   `p_meaningful_hi`, `p_equivalent`, `p_equivalent_mcse`, `p_equivalent_lo`,
#'   `p_equivalent_hi`, `n_attempted`, `n_returned`, `n_converged`, `n_singular`,
#'   `n_warning` and `fitter`. Each decision proportion is reported with its Monte Carlo standard
#'   error (`*_mcse`) and Wilson interval bounds (`*_lo`, `*_hi`), because a proportion
#'   over a finite number of replicates is an estimate rather than a fact. The interval
#'   behind `mean_ci_width` and the ROPE decisions is the Wald approximation described in
#'   Details. `fitter` names the fitter and the interval, such as
#'   `"lme4::glmer (binomial), Wald z"`. For an `ordinal` or `beta` response fitted on its own
#'   scale, the decision columns are `NA` (see Details).
#'
#'   The decision proportions are taken over `n_returned`, the replicates that
#'   produced an estimate. The remaining counts separate the fit outcomes, because a
#'   fit can return a usable estimate while still being boundary-singular or
#'   carrying a convergence warning: `n_converged` counts replicates with neither,
#'   `n_singular` those where `lme4::isSingular()` was true, and `n_warning` those
#'   with a warning or optimiser convergence message other than lme4's singular-fit
#'   notice, which `n_singular` already counts. Singular and warning fits are
#'   retained, since their fixed-effect estimates remain interpretable and
#'   discarding them would bias the result. A large `n_singular` means the model
#'   being fitted is richer than the design can support at that sample size, which
#'   is common in crossed designs (Bates et al., 2015; Matuschek et al., 2017).
#'
#'   A ladder adds one column per rung, `rung_<name>`, holding the share of the replicates
#'   analysed at that rung, the same on every row. The shares sum to 1, and a replicate in which no
#'   rung's fit met the rule counts at the last rung. The fit counts then describe the fit each
#'   replicate kept.
#' @references Bates, D., Kliegl, R., Vasishth, S. and Baayen, H. (2015). Parsimonious mixed
#'   models. \emph{arXiv}. \doi{10.48550/arXiv.1506.04967}
#'
#'   Matuschek, H., Kliegl, R., Vasishth, S., Baayen, H. and Bates, D. (2017). Balancing
#'   Type I error and power in linear mixed models. \emph{Journal of Memory and Language},
#'   94, 305-315. \doi{10.1016/j.jml.2017.01.001}
#'
#'   Bolker, B. M., Brooks, M. E., Clark, C. J., Geange, S. W., Poulsen, J. R., Stevens, M. H. H.
#'   and White, J.-S. S. (2009). Generalized linear mixed models: A practical guide for ecology
#'   and evolution. \emph{Trends in Ecology & Evolution}, 24(3), 127-135.
#'   \doi{10.1016/j.tree.2008.10.008}
#'
#'   Li, P. and Redden, D. T. (2015). Comparing denominator degrees of freedom approximations for
#'   the generalized linear mixed model in analyzing binary outcome in small sample
#'   cluster-randomized trials. \emph{BMC Medical Research Methodology}, 15, 38.
#'   \doi{10.1186/s12874-015-0026-x}
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
                             n_sims = 100, workers = 1, cl = NULL,
                             fallback = c("singular", "warning")) {
  if (!requireNamespace("lme4", quietly = TRUE))
    stop("precision_design() requires the 'lme4' package; please install it.", call. = FALSE)
  fallback <- match.arg(fallback)
  .refuse_focal_sets(focal, "precision_design()")
  spec <- .as_spec(spec)
  workers <- .check_workers(workers)
  cl <- .check_cluster(cl)
  rope <- .as_rope(rope)
  # A ladder is read and checked here, once, and travels to the replicate loop as the formula.
  ladder <- .as_ladder(formula, fallback, spec)
  if (!is.null(ladder)) formula <- ladder
  if (is.null(cl) && workers > 1L) {
    cl <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }
  .precision_design_impl(spec, focal, formula = formula, prep = prep, rope = rope,
                         n_sims = n_sims, cl = cl)
}

# The replicate loop behind precision_design(), run on the PSOCK cluster `cl` when there is one,
# which precision_curve() hands to every grid point through sweep_spec(). `rope` arrives as the
# region c(lower, upper) that .as_rope() returns. `formula` is NULL, a single formula or a ladder.
.precision_design_impl <- function(spec, focal, formula, prep, rope, n_sims, cl = NULL) {
  ladder <- if (inherits(formula, "pilotr_ladder")) formula
  # The response family chooses the fitter for the model pilotr derives, a ladder's structures
  # included, and a formula the user gives is fitted as written.
  auto <- is.null(formula) || isTRUE(ladder$derived)
  if (is.null(formula)) formula <- model_formula(spec)
  if (is.null(prep)) prep <- .default_prep(spec)
  fo <- .resolve_focal(focal, spec)
  fnames <- fo$names
  if (!length(fnames))
    stop("this specification has no fixed coefficients to analyse; name the focal effects explicitly",
         call. = FALSE)
  family <- spec[["response"]][["family"]]
  # An ordinal or beta response is fitted on its own scale, where the region of practical
  # equivalence, set on the coefficients' logit scale, means nothing, so the decisions are
  # withheld. The interval width stays, on the scale of the fitted model.
  withheld <- .linear_on_link(family, auto)
  if (withheld)
    warning(sprintf(.LINEAR_PRECISION, family, .LINK_SCALE[[family]]), call. = FALSE)

  seeds <- replicate_seeds(spec[["seed"]], n_sims)
  # The shared replicate loop, with the cheaper fitter: a precision analysis needs estimates and
  # standard errors, and lmerTest's Satterthwaite p-values cost noticeably more than the plain fit.
  res <- .p_lapply(seq_len(n_sims), .design_rep, cl = cl, spec = spec, seeds = seeds,
                   prep = prep, formula = formula, fnames = fnames, test = FALSE, auto = auto)

  n_returned <- sum(vapply(res, function(r) isTRUE(r$fitted), logical(1)))
  coef_names <- NULL
  for (r in res) if (isTRUE(r$fitted)) { coef_names <- r$coef_names; break }

  seen <- stats::setNames(integer(length(fnames)), fnames)
  blank <- stats::setNames(rep(NA_real_, length(fnames)), fnames)
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
    if (withheld) next
    d <- .rope_decisions(lo, hi, rope)
    ro <- .rate_with_error(sum(d$meaningful), length(ok), "p")
    ri <- .rate_with_error(sum(d$equivalent), length(ok), "p")
    p_out[f] <- ro$p; out_mcse[f] <- ro$p_mcse; out_lo[f] <- ro$p_lo; out_hi[f] <- ro$p_hi
    p_ins[f] <- ri$p; ins_mcse[f] <- ri$p_mcse; ins_lo[f] <- ri$p_lo; ins_hi[f] <- ri$p_hi
  }
  .warn_no_fits(res, n_returned, .last_rung(formula))
  .warn_absent_focal(seen, n_returned, coef_names)

  cnt <- .fit_counts(res, n_sims, n_returned)
  out <- data.frame(
    param = fnames, true = unname(fo$true[fnames]),
    mean_ci_width = unname(width),
    p_meaningful = unname(p_out), p_meaningful_mcse = unname(out_mcse),
    p_meaningful_lo = unname(out_lo), p_meaningful_hi = unname(out_hi),
    p_equivalent = unname(p_ins), p_equivalent_mcse = unname(ins_mcse),
    p_equivalent_lo = unname(ins_lo), p_equivalent_hi = unname(ins_hi),
    n_attempted = cnt$n_attempted, n_returned = cnt$n_returned,
    n_converged = cnt$n_converged, n_singular = cnt$n_singular, n_warning = cnt$n_warning,
    fitter = .analysis_label(formula, family, test = FALSE, auto = auto),
    row.names = NULL, stringsAsFactors = FALSE)
  # A ladder adds the share of the replicates at each rung, and a single model adds nothing.
  .add_rung_columns(out, if (!is.null(ladder)) .rung_shares(res, ladder))
}

# The default data-preparation function, built in its own small frame so that only the
# spec (not the caller's whole frame, cluster included) travels to PSOCK workers.
.default_prep <- function(spec) function(d) model_data(spec, d)

# The region of practical equivalence as c(lower, upper). A single finite half-width r >= 0 is
# the region c(-r, r), and a pair of finite numbers in increasing order is the region it names.
# Anything else is refused. hpc/precision_array.R hands its ROPE variable to precision_design(),
# which checks it here.
.as_rope <- function(rope) {
  ok <- is.numeric(rope) && all(is.finite(rope)) &&
    ((length(rope) == 1L && rope >= 0) || (length(rope) == 2L && rope[1] < rope[2]))
  if (!ok) {
    # Shown as code, so that NA, a string and a vector of the wrong length each read as what they
    # are. A long vector is cut short.
    shown <- paste(deparse(rope, width.cutoff = 500L), collapse = " ")
    if (nchar(shown) > 60L) shown <- paste0(substr(shown, 1L, 57L), "...")
    stop("`rope` must be a non-negative half-width or an increasing pair c(lower, upper); got ",
         shown, call. = FALSE)
  }
  rope <- as.numeric(rope)
  if (length(rope) == 1L) c(-rope, rope) else rope
}

# Each interval's decisions against the region c(lower, upper) from .as_rope(): meaningful when
# the interval lies wholly above or below the region, equivalent when it lies wholly inside. Both
# inequalities are strict, so an interval that touches a bound counts as neither. For a half-width
# r, they are the comparisons with -r and r. hpc/precision_array.R calls precision_design(), so it
# makes the same decisions.
.rope_decisions <- function(lo, hi, rope) {
  list(meaningful = lo > rope[2] | hi < rope[1],
       equivalent = lo > rope[1] & hi < rope[2])
}

#' Precision and ROPE curve over sample size
#'
#' Sweep the number of subjects and report the ROPE decision probabilities at each size. Pass the
#' result to [solve_curve()] for the minimum analysable \emph{N} at which a focal effect reaches a
#' determinate decision with a target probability, such as 0.90, together with an interval on it.
#' Calls [precision_design()] and so requires the `lme4` package.
#'
#' @details
#' A thin wrapper around [sweep_spec()] over `units$subject$n`, kept because a sample-size curve is
#' the sweep users want most often. For any other axis, call [sweep_spec()] directly; effect size is
#' the axis a design analysis most often needs next, and [design_conditions()] builds the
#' coefficient overrides for it.
#'
#' Each sample size draws its own replicates, from the seed [sweep_spec()] gives its grid point,
#' so the points of the curve are independent, as [solve_curve()] assumes when it fits them.
#' `seeds = "shared"` runs every point from the specification's own seed instead, and the points
#' then share their noise (see [sweep_spec()]).
#' With `workers > 1`, one worker pool is created and reused across all sample sizes.
#'
#' Reading the crossing off the returned curve, or off a plot of it, is what [solve_curve()]
#' replaces. At the replicate counts these runs are usually given, neighbouring points on the
#' curve are not significantly different from one another, so a sample size judged by eye is a
#' point estimate with an unstated and often wide uncertainty behind it.
#'
#' @param spec A design specification (path or list).
#' @param focal The focal effects, as in [precision_design()]. `NULL` uses every coefficient in the
#'   specification.
#' @param subject_ns A numeric vector of subject counts to evaluate.
#' @param formula Optional `lme4` formula; if `NULL` it is derived via
#'   [model_formula()], and the response family chooses the fitter, as in [precision_design()]. A
#'   ladder of models, as in [precision_design()], simulates a simplification rule at every sample
#'   size.
#' @param prep Optional data-preparation function; if `NULL` it is derived via
#'   [model_data()].
#' @param rope The region of practical equivalence, as in [precision_design()]: a non-negative
#'   half-width `r` for the region from `-r` to `r`, or a pair `c(lower, upper)` with
#'   `lower < upper`. Keep the effects worth detecting clearly outside the region, because the
#'   probability of a determinate meaningful decision about an effect inside it cannot rise above
#'   0.5 however large the sample, so the curve would fall with \emph{N} rather than rise.
#' @param n_sims Number of Monte Carlo replicates per sample size. `p_meaningful` and
#'   `p_equivalent` are proportions over the replicates that produced an estimate, so each is
#'   reported with its Monte Carlo standard error and Wilson interval. The default of 60 gives a
#'   standard error of 0.065 at a rate of 0.5, which is too coarse to support a claim about a
#'   design; raise it to at least 200 for planning.
#' @param workers Number of local worker processes over which to spread the replicates at
#'   each sample size. The default of 1 runs serially, and any worker count returns results
#'   identical to a serial run. One pool of that many processes serves the whole sweep.
#' @param seeds How the sample sizes are seeded: `"independent"`, the default, or `"shared"`, as
#'   in [sweep_spec()]'s `.seeds`.
#' @param fallback The rule of a ladder, `"singular"` or `"warning"`, as in [precision_design()].
#' @return A data frame with one row per focal effect and sample size, adding an `n_subject`
#'   column to the columns returned by [precision_design()], including the Monte Carlo standard
#'   errors, the Wilson interval bounds, and the `n_returned`, `n_converged`, `n_singular` and
#'   `n_warning` fit counts, and for a ladder the share of each rung at each sample size.
#' @seealso [solve_curve()] to solve the returned curve for a target decision probability,
#'   [sweep_spec()] for any other axis, and [design_conditions()] for effect-size grids.
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
                            rope = 0.05, n_sims = 60, workers = 1,
                            seeds = c("independent", "shared"),
                            fallback = c("singular", "warning")) {
  if (!requireNamespace("lme4", quietly = TRUE))
    stop("precision_curve() requires the 'lme4' package; please install it.", call. = FALSE)
  .refuse_focal_sets(focal, "precision_curve()")
  spec <- .as_spec(spec)
  seeds <- match.arg(seeds)
  fallback <- match.arg(fallback)
  workers <- .check_workers(workers)
  .as_rope(rope)                    # a malformed region is refused before any worker starts
  .as_ladder(formula, fallback, spec)   # and so is a malformed ladder
  # One pool for the whole sweep, handed to every grid point, so that precision_design() does
  # not start and stop a pool at every point.
  cl <- NULL
  if (workers > 1L) {
    cl <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }
  sweep_spec(spec, "units$subject$n", subject_ns, precision_design,
             focal = focal, formula = formula, prep = prep, rope = rope,
             n_sims = n_sims, cl = cl, fallback = fallback, .name = "n_subject", .seeds = seeds)
}
