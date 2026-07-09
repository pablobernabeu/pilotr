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
#' @param focal A named numeric vector mapping focal coefficient names to their true values,
#'   or a character vector of coefficient names.
#' @param formula Optional `lme4` formula; if `NULL` it is derived from the specification via
#'   [model_formula()].
#' @param prep Optional function mapping a simulated data set to the modelling data; if `NULL`
#'   it is derived via [model_data()], which log-transforms the outcome and builds the
#'   contrast and interaction columns, so focal names follow the auto-formula (interactions
#'   written as `a_b`).
#' @param rope Half-width of the region of practical equivalence; an effect with
#'   `abs(beta) < rope` is treated as practically equivalent to zero.
#' @param n_sims Number of Monte Carlo replicates.
#' @param workers Number of local worker processes over which to spread the replicates.
#'   The default of 1 runs serially. Because every replicate seeds the shared RNG from its
#'   own index, any worker count returns results identical to a serial run.
#' @return A data frame with one row per focal effect and columns `param`, `true`,
#'   `mean_ci_width`, `p_meaningful`, `p_equivalent`, and `n_converged`. The interval
#'   behind `mean_ci_width` and the ROPE decisions is the Wald approximation described in
#'   Details.
#' @examples
#' \donttest{
#' if (requireNamespace("lme4", quietly = TRUE)) {
#'   spec <- build_spec(list(name = "pr", seed = 1, design_kind = "within",
#'     include_items = TRUE, n_subject = 12, n_item = 12, factor_name = "cond",
#'     lev1 = "a", lev2 = "b", intercept = 6, effect = 0.05,
#'     subj_int_sd = 0.12, subj_slope_sd = 0.04, subj_corr = 0.2,
#'     item_int_sd = 0.08, item_slope_sd = 0.02, item_corr = -0.1,
#'     family = "shifted_lognormal", resp_name = "", sigma = 0.3, shift = 200))
#'   precision_design(spec, focal = c(effect = 0.05), rope = 0.02, n_sims = 10)
#' }
#' }
#' @export
precision_design <- function(spec, focal, formula = NULL, prep = NULL, rope = 0.05,
                             n_sims = 100, workers = 1) {
  if (!requireNamespace("lme4", quietly = TRUE))
    stop("precision_design() requires the 'lme4' package; please install it.", call. = FALSE)
  if (is.character(spec)) spec <- load_spec(spec)
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
  fnames <- if (!is.null(names(focal))) names(focal) else as.character(focal)
  W <- setNames(lapply(fnames, function(.) numeric(0)), fnames)
  OUT <- setNames(integer(length(fnames)), fnames); INS <- OUT; NC <- 0L
  base <- spec$seed
  res <- .p_lapply(seq_len(n_sims), .precision_rep, cl = cl, spec = spec, base = base,
                   prep = prep, formula = formula, fnames = fnames, rope = rope)
  for (r in res) {
    if (is.null(r)) next
    NC <- NC + 1L
    for (f in fnames) {
      if (!r$present[[f]]) next
      W[[f]] <- c(W[[f]], r$width[[f]])
      if (r$out[[f]]) OUT[f] <- OUT[f] + 1L
      if (r$ins[[f]]) INS[f] <- INS[f] + 1L
    }
  }
  data.frame(
    param = fnames,
    true = if (!is.null(names(focal))) unname(focal) else NA_real_,
    mean_ci_width = vapply(fnames, function(f) if (length(W[[f]])) mean(W[[f]]) else NA_real_, numeric(1)),
    p_meaningful = OUT / NC, p_equivalent = INS / NC, n_converged = NC, row.names = NULL)
}

# The default data-preparation function, built in its own small frame so that only the
# spec (not the caller's whole frame, cluster included) travels to PSOCK workers.
.default_prep <- function(spec) function(d) model_data(spec, d)

# One precision replicate. Kept at top level so that only the arguments travel to PSOCK
# workers. Returns NULL when the fit fails; otherwise, per focal effect, the CI width and
# the two ROPE decisions, with `present` recording whether the effect was in the fit.
.precision_rep <- function(i, spec, base, prep, formula, fnames, rope) {
  s <- spec; s$seed <- base + (i - 1)          # same indexed-seed rule as the power functions
  d <- prep(simulate_design(s))
  fit <- tryCatch(suppressWarnings(suppressMessages(
    lme4::lmer(formula, data = d, control = lme4::lmerControl(calc.derivs = FALSE)))),
    error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  est <- lme4::fixef(fit); se <- sqrt(diag(as.matrix(stats::vcov(fit))))
  present <- setNames(fnames %in% names(est), fnames)
  width <- setNames(rep(NA_real_, length(fnames)), fnames)
  out <- setNames(logical(length(fnames)), fnames); ins <- out
  for (f in fnames[present]) {
    lo <- est[f] - 1.96 * se[f]; hi <- est[f] + 1.96 * se[f]
    width[f] <- hi - lo
    out[f] <- lo > rope || hi < -rope
    ins[f] <- lo > -rope && hi < rope
  }
  list(present = present, width = width, out = out, ins = ins)
}

#' Precision and ROPE curve over sample size
#'
#' Sweep the number of subjects and report the ROPE decision probabilities at each size, to
#' identify the minimum analysable \emph{N} at which a focal effect reaches a determinate decision
#' with a target probability (for example 0.90). Calls [precision_design()] and so requires
#' the `lme4` package.
#'
#' @param spec A design specification (path or list).
#' @param focal A named numeric vector (coefficient name to true value) or character vector
#'   of focal coefficient names.
#' @param subject_ns A numeric vector of subject counts to evaluate.
#' @param formula Optional `lme4` formula; if `NULL` it is derived via [model_formula()].
#' @param prep Optional data-preparation function; if `NULL` it is derived via [model_data()].
#' @param rope Half-width of the region of practical equivalence.
#' @param n_sims Number of Monte Carlo replicates per sample size.
#' @param workers Number of local worker processes over which to spread the replicates at
#'   each sample size. The default of 1 runs serially, and any worker count returns results
#'   identical to a serial run. The worker pool is created once and reused across the sweep.
#' @return A data frame with one row per focal effect and sample size, adding an `n_subject`
#'   column to the columns returned by [precision_design()].
#' @examples
#' \donttest{
#' if (requireNamespace("lme4", quietly = TRUE)) {
#'   spec <- build_spec(list(name = "pr", seed = 1, design_kind = "within",
#'     include_items = TRUE, n_subject = 12, n_item = 12, factor_name = "cond",
#'     lev1 = "a", lev2 = "b", intercept = 6, effect = 0.05,
#'     subj_int_sd = 0.12, subj_slope_sd = 0.04, subj_corr = 0.2,
#'     item_int_sd = 0.08, item_slope_sd = 0.02, item_corr = -0.1,
#'     family = "shifted_lognormal", resp_name = "", sigma = 0.3, shift = 200))
#'   precision_curve(spec, focal = c(effect = 0.05), subject_ns = c(12, 18), n_sims = 8)
#' }
#' }
#' @export
precision_curve <- function(spec, focal, subject_ns, formula = NULL, prep = NULL, rope = 0.05,
                            n_sims = 60, workers = 1) {
  if (!requireNamespace("lme4", quietly = TRUE))
    stop("precision_curve() requires the 'lme4' package; please install it.", call. = FALSE)
  if (is.character(spec)) spec <- load_spec(spec)
  if (is.null(formula)) formula <- model_formula(spec)
  if (is.null(prep)) prep <- .default_prep(spec)
  workers <- .check_workers(workers)
  cl <- NULL
  if (workers > 1L) {
    cl <- parallel::makeCluster(workers)   # one cluster, reused across the whole sweep
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }
  parts <- lapply(subject_ns, function(n) {
    s <- spec; s$units$subject$n <- n
    cbind(n_subject = n, .precision_design_impl(s, focal, formula = formula, prep = prep,
                                                rope = rope, n_sims = n_sims, cl = cl))
  })
  do.call(rbind, parts)
}
