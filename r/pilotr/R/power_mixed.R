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
#' specification implies and test each focal fixed effect. The fitter follows the response family:
#' `lmerTest::lmer()` with Satterthwaite p-values for the continuous families, and `lme4::glmer()`
#' with a Wald z for accuracy and count responses (see Details). Reports power together with the
#' Type S and Type M errors of Gelman and Carlin (2014). Requires the `lme4` and `lmerTest`
#' packages.
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
#' The fitter follows the response family, so that the estimates are on the scale of the
#' specification's coefficients. A `gaussian`, `lognormal`, `shifted_lognormal` or `exgaussian`
#' response is fitted by `lmerTest::lmer()`, the two lognormal families on the log scale to which
#' [model_data()] takes them, and each effect is tested with Satterthwaite's degrees of freedom. A
#' `bernoulli` or `poisson` response is fitted by `lme4::glmer()` with the binomial or Poisson
#' family, on the logit or log scale, and each effect is tested with a Wald z. A model with no
#' random terms, such as that of a between-subjects design with one row per subject, is fitted by
#' `stats::lm()` or `stats::glm()` in the same way. The result's `fitter` names the fitter and the
#' test. A Wald z treats the estimate over its standard error as normal and so ignores the
#' uncertainty in the variance components. Bolker et al. (2009) recommend it for the fixed effects
#' of a GLMM without overdispersion, which these models are, but it rejects too often when the
#' clusters are few. In binary GLMMs with fewer than 30 clusters, a test with as many degrees of
#' freedom as observations, close to a z test, rejected more often than its nominal level (Li and
#' Redden, 2015). The power of a `bernoulli` or `poisson` design with few subjects or items is
#' therefore somewhat overstated.
#'
#' pilotr has no frequentist model for an `ordinal` or `beta` response yet. Such a response is
#' fitted by the linear model on its own scale, while the specification's coefficients are on the
#' logit scale. The function therefore warns, and withholds the mean estimate and Type M (`NA`),
#' since both would compare the two scales. Type S compares signs alone, and the link keeps the
#' sign of a main effect unless an interaction reverses that effect between the levels of another
#' factor. It does not keep the sign of an interaction. A bounded response compresses differences
#' near the ends of its range, so an interaction that is positive on the logit scale can be
#' negative on the response scale. Type S is therefore withheld for an interaction as well. Power
#' is the rate at which the linear model's test rejects. [generate_design_analysis()] writes a
#' Bayesian design analysis on the link scale, with the model [brms_bridge()] derives.
#'
#' A `formula` given here is fitted as written, by `lmerTest::lmer()`, or by `stats::lm()` when
#' it has no random terms, whatever the family, since only its author knows the scale of its
#' response. With such a formula, or a `prep` that changes the response, and `focal = NULL`, the
#' true values still come from the specification and are on its scale, which need not be the
#' model's. A named numeric `focal` gives the true values on the model's scale.
#'
#' Every reported rate carries its Monte Carlo standard error and a Wilson interval, because a
#' proportion over a finite number of replicates is an estimate rather than a fact. At the default
#' 100 replicates a power near 0.5 has a standard error of 0.05. Type S and Type M average over
#' the significant replicates alone, which may be few. Each carries its own Monte Carlo error.
#' Type S has a Wilson interval, since its standard error would be zero whenever no sign is wrong.
#' Type M has its standard error, the standard deviation of the exaggeration ratios over the square
#' root of their number.
#'
#' A directional prediction calls for a one-sided test, which `alternative` selects for every
#' focal effect. `"greater"` tests whether the coefficient is positive and `"less"` whether it is
#' negative. The one-sided p-value is the upper or lower tail of the statistic behind the
#' two-sided one. That is Satterthwaite's t, the t of `stats::lm()` with its residual degrees of
#' freedom, or a Wald z referred to the normal distribution. The one-sided p-value is half the
#' two-sided one when the estimate lies in the tested direction. A one-sided test at `alpha`
#' therefore rejects in exactly the replicates in which the two-sided test at `2 * alpha` rejects
#' with the estimate in that direction. Type S and Type M are taken over the replicates
#' significant in the tested direction. When that is the direction of the true effect, Type S is
#' 0 by construction, and Type M is the exaggeration among the replicates that would support the
#' prediction. `fitter` names the direction, as in
#' `"lmerTest::lmer, Satterthwaite t, one-sided (greater)"`.
#'
#' Results from earlier versions were two-sided. A one-sided power at `alpha`, in the direction of
#' the true effect, can be read off a two-sided result at `2 * alpha`: it is
#' `power * (1 - type_s)`. Type M is recovered only when that result's Type S is 0, since the
#' wrong-signed replicates cannot otherwise be told apart from the others.
#'
#' @param spec A design specification (path or list).
#' @param focal The fixed effects to test. `NULL`, the default, tests every coefficient in the
#'   specification and takes the true values from it. A character vector names the effects and
#'   leaves the true values unknown, which suppresses Type S and Type M. A named numeric vector
#'   gives both, which is how to test against a value other than the one simulated. Interaction
#'   effects follow the model's column naming, so a specification key `a:b` is the focal name
#'   `a_b`.
#' @param formula Optional `lme4` formula; if `NULL` it is derived from the specification via
#'   [model_formula()], and the response family chooses the fitter. A formula given here is
#'   fitted as a linear model whatever the family (see Details).
#' @param prep Optional function mapping a simulated data set to the modelling data; if `NULL` it
#'   is derived via [model_data()].
#' @param n_sims Number of Monte Carlo replicates. A power estimate carries a Monte Carlo
#'   standard error of about `sqrt(p * (1 - p) / n_sims)`, and `type_s` and
#'   `type_m` average over the significant replicates alone, so they settle more
#'   slowly still. At least 200 replicates are advisable for study planning.
#' @param alpha Significance level of the test, two-sided unless `alternative` says otherwise.
#' @param workers Number of local worker processes over which to spread the replicates.
#'   The default of 1 runs serially. Because the replicate seeds are derived once from the
#'   specification's seed, any worker count returns results identical to a serial run. The
#'   mixed-model fits dominate the cost, so the speed-up is close to linear in the number of cores.
#' @param cl A cluster made by [parallel::makeCluster()] to run the replicates on, in place of one
#'   started for this call, in which case `workers` is ignored. The cluster is left running for its
#'   owner to stop. Through [sweep_spec()], as in
#'   `sweep_spec(spec, "fixed$coefficients$cond", values, power_mixed, cl = cl)`, one cluster serves
#'   every grid point.
#' @param alternative The alternative hypothesis of every focal effect's test: `"two.sided"`, the
#'   default, or `"greater"` or `"less"` for a one-sided test of a directional prediction (see
#'   Details).
#' @return An object of class `pilotr_power`, a list whose per-run elements are `n_sims`,
#'   `alpha`, `fitter`, `n_attempted`, `n_returned`, `n_converged`, `n_singular` and `n_warning`,
#'   and whose per-effect elements are vectors named by focal effect: `power`, `power_mcse`,
#'   `power_lo`, `power_hi`, `n_significant`, `true_effect`, `mean_estimate`, `type_s`,
#'   `type_s_lo`, `type_s_hi`, `type_m` and `type_m_mcse`. With a single focal effect each of
#'   those has length one, so `result$power` reads as it always has. `fitter` names the fitter and
#'   its test, such as `"lme4::glmer (binomial), Wald z"`, with the direction of a one-sided test
#'   added, as in `"lme4::glmer (binomial), Wald z, one-sided (greater)"`. Under a one-sided
#'   test, `n_significant` and the design-analysis quantities count the replicates significant in
#'   the tested direction. `mean_estimate` and `type_m` are `NA` for an `ordinal` or `beta`
#'   response fitted on its own scale, and so is the `type_s` of an interaction (see Details).
#'   `type_s_lo` and `type_s_hi`, the Wilson interval of `type_s`, and `type_m_mcse`, the Monte
#'   Carlo standard error of `type_m`, are `NA` wherever their rate is, and `type_m_mcse` also
#'   when a single replicate reached significance.
#'
#'   `power` is the proportion of significant results among the replicates that returned an
#'   estimate for that effect, not among `n_sims`. The counts report the fit outcomes separately,
#'   because a fit can return a usable estimate while still being boundary-singular or carrying a
#'   convergence warning: `n_returned` counts replicates that yielded a fit, `n_converged` those
#'   that did so with neither a warning nor a singular fit, `n_singular` those where
#'   `lme4::isSingular()` was true, and `n_warning` those with a warning or optimiser convergence
#'   message other than lme4's singular-fit notice, which `n_singular` already counts. Singular and warning fits are retained in `power`, since their fixed-effect estimates
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
                        n_sims = 100, alpha = 0.05, workers = 1, cl = NULL,
                        alternative = c("two.sided", "greater", "less")) {
  if (!requireNamespace("lme4", quietly = TRUE) ||
      !requireNamespace("lmerTest", quietly = TRUE))
    stop("power_mixed() requires the 'lme4' and 'lmerTest' packages; please install them.",
         call. = FALSE)
  alternative <- match.arg(alternative)
  spec <- .as_spec(spec)
  workers <- .check_workers(workers)
  cl <- .check_cluster(cl)
  if (is.null(cl) && workers > 1L) {
    cl <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }
  .power_mixed_impl(spec, focal = focal, formula = formula, prep = prep,
                    n_sims = n_sims, alpha = alpha, cl = cl, alternative = alternative)
}

# The replicate loop behind power_mixed(), run on the PSOCK cluster `cl` when there is one. The
# cluster is the caller's or one power_mixed() started, and power_curve_mixed() hands one to every
# grid point through sweep_spec().
.power_mixed_impl <- function(spec, focal, formula, prep, n_sims, alpha, cl = NULL,
                              alternative = "two.sided") {
  auto <- is.null(formula)
  if (auto) formula <- model_formula(spec)
  if (is.null(prep)) prep <- .default_prep(spec)
  fo <- .resolve_focal(focal, spec)
  fnames <- fo$names
  if (!length(fnames))
    stop("this specification has no fixed coefficients to test; name the focal effects explicitly",
         call. = FALSE)
  family <- spec[["response"]][["family"]]
  # An ordinal or beta response is fitted on its own scale while the true values are on the logit
  # scale, so the comparisons between the two are withheld. The warning is raised here, once per
  # call, since the replicates may run in other processes.
  withheld <- .linear_on_link(family, auto)
  if (withheld)
    warning(sprintf(.LINEAR_POWER, family, .LINK_SCALE[[family]]), call. = FALSE)
  # The link keeps the sign of a difference between two cells, but not that of an interaction, a
  # difference of differences. A bounded response compresses differences near the ends of its
  # range, so a positive interaction on the logit scale can be negative on the response scale.
  # On a beta response with cell means of 0.525, 0.802, 0.802 and 0.957, a logit-scale interaction
  # of 0.4 is -0.122 on the response scale, and Type S read 1 in every replicate.
  interactions <- if (withheld) .us(.interaction_keys(spec)) else character(0)

  seeds <- replicate_seeds(spec[["seed"]], n_sims)
  res <- .p_lapply(seq_len(n_sims), .design_rep, cl = cl, spec = spec, seeds = seeds,
                   prep = prep, formula = formula, fnames = fnames, test = TRUE, auto = auto,
                   alternative = alternative)

  n_returned <- sum(vapply(res, function(r) isTRUE(r$fitted), logical(1)))
  coef_names <- NULL
  for (r in res) if (isTRUE(r$fitted)) { coef_names <- r$coef_names; break }

  seen <- stats::setNames(integer(length(fnames)), fnames)
  blank <- stats::setNames(rep(NA_real_, length(fnames)), fnames)
  power <- mcse <- lo <- hi <- type_s <- s_lo <- s_hi <- type_m <- m_mcse <- mean_est <- blank
  n_sig <- stats::setNames(integer(length(fnames)), fnames)

  for (f in fnames) {
    est <- vapply(res, function(r) r$est[[f]], numeric(1))
    pv <- vapply(res, function(r) r$p[[f]], numeric(1))
    ok <- which(!is.na(pv))
    seen[f] <- length(ok)
    if (!length(ok)) next
    # Under a one-sided test, the replicates significant in the tested direction, so Type S and
    # Type M below are taken over them alone.
    sig <- ok[pv[ok] < alpha]
    rate <- .rate_with_error(length(sig), length(ok), "power")
    power[f] <- rate$power; mcse[f] <- rate$power_mcse
    lo[f] <- rate$power_lo; hi[f] <- rate$power_hi
    n_sig[f] <- length(sig)
    if (!withheld) mean_est[f] <- mean(est[ok])
    beta <- fo$true[[f]]
    # Type S and Type M are defined relative to a true value, and Type M divides by it, so both
    # stay NA when the true effect is unknown or zero rather than reporting an infinity. Type S
    # compares signs alone, which a monotone link keeps for a main effect, so it survives a fit on
    # the response scale except for an interaction. Type M compares magnitudes, which the link
    # changes. Each carries its Monte Carlo error and is withheld together with it: a Wilson
    # interval for Type S and a standard error for Type M, which needs two significant replicates.
    if (length(sig) && !is.na(beta) && beta != 0) {
      if (!(f %in% interactions)) {
        wrong <- (est[sig] > 0) != (beta > 0)
        type_s[f] <- mean(wrong)
        ci <- .type_s_interval(sum(wrong), length(sig), beta)
        s_lo[f] <- ci[1]; s_hi[f] <- ci[2]
      }
      if (!withheld) {
        ratios <- abs(est[sig]) / abs(beta)
        type_m[f] <- mean(ratios)
        m_mcse[f] <- .type_m_mcse(ratios)
      }
    }
  }
  .warn_no_fits(res, n_returned, formula)
  .warn_absent_focal(seen, n_returned, coef_names)

  out <- c(list(n_sims = n_sims, alpha = alpha,
                fitter = .fitter_label(formula, family, test = TRUE, auto = auto,
                                       alternative = alternative)),
           .fit_counts(res, n_sims, n_returned),
           list(power = power, power_mcse = mcse, power_lo = lo, power_hi = hi,
                n_significant = n_sig, true_effect = fo$true[fnames],
                mean_estimate = mean_est, type_s = type_s, type_s_lo = s_lo, type_s_hi = s_hi,
                type_m = type_m, type_m_mcse = m_mcse))
  class(out) <- c("pilotr_power", "list")
  out
}

#' Print a simulation-based power result
#'
#' Shows the fit accounting and then two tables with one row per focal effect. The first gives
#' each power estimate beside its Monte Carlo standard error and Wilson interval, so that the
#' precision of the estimate is as visible as the estimate. The second gives Type S beside its
#' Wilson interval and Type M beside its Monte Carlo standard error, both taken over the
#' significant replicates alone.
#'
#' @param x A `pilotr_power` object, as returned by [power_mixed()].
#' @param digits Number of decimal places for the reported rates and for the interval bounds.
#' @param ... Ignored, present for consistency with the generic.
#' @return `x`, invisibly.
#' @export
print.pilotr_power <- function(x, digits = 3, ...) {
  cat(sprintf("Simulation-based power over %d replicates (alpha = %g)\n", x$n_sims, x$alpha))
  # A result saved by a version that did not record its fitter prints as it did.
  if (!is.null(x$fitter)) cat(sprintf("  fitter: %s\n", x$fitter))
  cat(sprintf("  fits: %d attempted, %d returned, %d converged cleanly, %d singular, %d with warnings\n",
              x$n_attempted, x$n_returned, x$n_converged, x$n_singular, x$n_warning))
  if (x$n_singular > 0.25 * max(x$n_attempted, 1))
    cat("  note: many fits were boundary-singular, so this model is richer than the design supports\n")
  # The bounds take `digits` as the rates do. They were written with three decimals whatever
  # `digits` said.
  interval <- function(lo, hi)
    ifelse(is.na(lo), NA_character_, sprintf("[%.*f, %.*f]", digits, lo, digits, hi))
  effect <- names(x$power)
  print(data.frame(
    effect = effect,
    true = unname(x$true_effect),
    power = round(unname(x$power), digits),
    mcse = round(unname(x$power_mcse), digits),
    ci95 = interval(unname(x$power_lo), unname(x$power_hi)),
    n_sig = unname(x$n_significant),
    row.names = NULL, stringsAsFactors = FALSE), row.names = FALSE)
  # Only the n_sig significant replicates inform Type S and Type M, so each is shown with its
  # uncertainty, in a second table that keeps both within 80 columns. A result saved by a version
  # that did not report the uncertainty prints the two rates alone.
  cat("  Type S and Type M, over the significant replicates:\n")
  design <- data.frame(effect = effect, type_s = round(unname(x$type_s), digits),
                       row.names = NULL, stringsAsFactors = FALSE)
  if (!is.null(x$type_s_lo))
    design$type_s_ci95 <- interval(unname(x$type_s_lo), unname(x$type_s_hi))
  design$type_m <- round(unname(x$type_m), digits)
  if (!is.null(x$type_m_mcse))
    design$type_m_mcse <- round(unname(x$type_m_mcse), digits)
  print(design, row.names = FALSE)
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
#' Each sample size draws its own replicates, from the seed [sweep_spec()] gives its grid point,
#' so the points of the curve are independent, as [target_n()] assumes when it fits them.
#' `seeds = "shared"` runs every point from the specification's own seed instead, as up to 0.3.1.
#'
#' For any axis other than sample size, call [sweep_spec()] directly. Effect size is the axis a
#' design analysis most often needs after sample size, and [design_conditions()] builds the
#' coefficient overrides for it.
#'
#' The curve is the input to [target_n()], which is where the sample size a preregistration
#' quotes should come from. Reading the crossing off a plot instead judges points whose Monte
#' Carlo intervals overlap, and reports the answer without the interval that goes with it.
#'
#' In a crossed design, the number of items limits what more subjects can buy. When the effect of
#' a condition varies from item to item, as a by-item random slope says it does, part of the
#' standard error of the effect comes from the items sampled. That part shrinks only as items are
#' added. As subjects are added, the standard error therefore falls towards a floor set by the
#' number of items. Power at a fixed number of items levels off below 1 (Westfall, Kenny and Judd,
#' 2014; Judd, Westfall and Kenny, 2017). The probit curve that [target_n()] fits rises towards 1
#' and cannot bend into such a plateau. A curve that levels off near the target can then be solved far
#' from the true crossing. One that levels off just below the target never crosses it, but Monte
#' Carlo noise can lift a simulated power over the target, and the curve can then be solved as
#' though it crossed. To compare designs that trade subjects for items, sweep `units$item$n` with
#' [sweep_spec()] at each of a few subject counts and solve each curve for the number of items, as
#' the power vignette shows. Power over items levels off in the same way, at a ceiling set by the
#' by-subject slope and the number of subjects. The subject counts must be large enough to put that
#' ceiling well above the target.
#'
#' @param spec A design specification (path or list).
#' @param subject_ns A numeric vector of subject counts to evaluate.
#' @param focal The fixed effects to test, as in [power_mixed()].
#' @param n_sims Number of Monte Carlo replicates per point. A power estimate carries a Monte
#'   Carlo standard error of about `sqrt(p * (1 - p) / n_sims)`, reported alongside it as
#'   `power_mcse`. The default of 60 gives a standard error of 0.065 at a power of 0.5, which is
#'   too coarse to support a claim about a design; raise it to at least 200 for planning.
#' @param alpha Significance level of the test, two-sided unless `alternative` says otherwise.
#' @param workers Number of local worker processes over which to spread the replicates at
#'   each grid point. The default of 1 runs serially, and any worker count returns results
#'   identical to a serial run. One pool of that many processes serves the whole sweep.
#' @param seeds How the sample sizes are seeded: `"independent"`, the default, or `"shared"`, as
#'   in [sweep_spec()]'s `.seeds`.
#' @param alternative The alternative hypothesis of every test, `"two.sided"`, `"greater"` or
#'   `"less"`, as in [power_mixed()].
#' @return A data frame with one row per sample size and focal effect, with columns `n_subject`,
#'   `effect`, `true`, `power`, `power_mcse`, `power_lo`, `power_hi`, `n_significant`, `type_s`,
#'   `type_s_lo`, `type_s_hi`, `type_m`, `type_m_mcse`, the `n_attempted`, `n_returned`,
#'   `n_converged`, `n_singular` and `n_warning` fit counts, and `fitter`, as in [power_mixed()].
#'   `fitter` names the direction of a one-sided test.
#'   `n_singular` typically falls as the sample size rises, so reading it down the sweep shows
#'   where the model becomes supportable.
#' @references Green, P. and MacLeod, C. J. (2016). SIMR: An R package for power analysis
#'   of generalized linear mixed models by simulation. \emph{Methods in Ecology and
#'   Evolution}, 7(4), 493-498. \doi{10.1111/2041-210x.12504}
#'
#'   Judd, C. M., Westfall, J. and Kenny, D. A. (2017). Experiments with more than one random
#'   factor: Designs, analytic models, and statistical power. \emph{Annual Review of Psychology},
#'   68, 601-625. \doi{10.1146/annurev-psych-122414-033702}
#'
#'   Kumle, L., Vo, M. L.-H. and Draschkow, D. (2021). Estimating power in (generalized)
#'   linear mixed models: An open introduction and tutorial in R. \emph{Behavior Research
#'   Methods}, 53, 2528-2543. \doi{10.3758/s13428-021-01546-0}
#'
#'   Westfall, J., Kenny, D. A. and Judd, C. M. (2014). Statistical power and optimal design in
#'   experiments in which samples of participants respond to samples of stimuli. \emph{Journal of
#'   Experimental Psychology: General}, 143(5), 2020-2045. \doi{10.1037/xge0000014}
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
                              workers = 1, seeds = c("independent", "shared"),
                              alternative = c("two.sided", "greater", "less")) {
  if (!requireNamespace("lme4", quietly = TRUE) ||
      !requireNamespace("lmerTest", quietly = TRUE))
    stop("power_curve_mixed() requires the 'lme4' and 'lmerTest' packages; please install them.",
         call. = FALSE)
  spec <- .as_spec(spec)
  seeds <- match.arg(seeds)
  alternative <- match.arg(alternative)
  workers <- .check_workers(workers)
  # One pool for the whole sweep, handed to every grid point. Since 0.3.0, when the sweep moved
  # into sweep_spec(), `workers` had travelled to power_mixed(), which started and stopped a pool
  # at every point.
  cl <- NULL
  if (workers > 1L) {
    cl <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }
  sweep_spec(spec, "units$subject$n", subject_ns, power_mixed,
             focal = focal, n_sims = n_sims, alpha = alpha, cl = cl, alternative = alternative,
             .name = "n_subject", .seeds = seeds)
}

# Flatten a pilotr_power object into one data-frame row per focal effect, for the sweeps.
.as_power_frame <- function(x) {
  data.frame(
    effect = names(x$power), true = unname(x$true_effect),
    power = unname(x$power), power_mcse = unname(x$power_mcse),
    power_lo = unname(x$power_lo), power_hi = unname(x$power_hi),
    n_significant = unname(x$n_significant),
    type_s = unname(x$type_s), type_s_lo = unname(x$type_s_lo), type_s_hi = unname(x$type_s_hi),
    type_m = unname(x$type_m), type_m_mcse = unname(x$type_m_mcse),
    n_attempted = x$n_attempted, n_returned = x$n_returned, n_converged = x$n_converged,
    n_singular = x$n_singular, n_warning = x$n_warning,
    fitter = if (is.null(x$fitter)) NA_character_ else x$fitter,
    row.names = NULL, stringsAsFactors = FALSE)
}
