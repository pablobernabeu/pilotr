# Precision-based design analysis against a Region of Practical Equivalence (ROPE), together
# with an N-sweep to find the minimum analysable sample size. This is a fast frequentist
# analogue of a Bayesian HDI-versus-ROPE design analysis. Across Monte Carlo replicates we
# record, for each focal fixed effect, the 95% CI and whether it falls determinately outside
# the ROPE (a practically meaningful effect) or entirely inside it (practical equivalence).

#' Precision/ROPE design analysis at a fixed design size.
#'
#' @param spec a design spec (path or list).
#' @param formula an `lmer` formula evaluated on the prepared data.
#' @param focal named numeric vector (coefficient name -> true value) or character vector.
#' @param rope half-width of the ROPE (|beta| < rope is "practically equivalent to zero").
#' @param n_sims number of Monte Carlo replicates.
#' @param formula optional `lmer` formula; if NULL it is auto-derived via [model_formula()].
#' @param prep optional function mapping simulated data to modelling data. If NULL it is
#'   auto-derived via [model_data()], which log-transforms the outcome and builds the contrast
#'   and interaction columns. Focal names then follow the auto-formula (interactions as `a_b`).
#' @export
precision_design <- function(spec, focal, formula = NULL, prep = NULL, rope = 0.05, n_sims = 100) {
  if (is.character(spec)) spec <- load_spec(spec)
  if (is.null(formula)) formula <- model_formula(spec)
  if (is.null(prep)) prep <- function(d) model_data(spec, d)
  fnames <- if (!is.null(names(focal))) names(focal) else as.character(focal)
  W <- setNames(lapply(fnames, function(.) numeric(0)), fnames)
  OUT <- setNames(integer(length(fnames)), fnames); INS <- OUT; NC <- 0L
  base <- spec$seed
  for (i in seq_len(n_sims)) {
    s <- spec; s$seed <- base + i
    d <- prep(simulate_design(s))
    fit <- tryCatch(suppressWarnings(suppressMessages(
      lme4::lmer(formula, data = d, control = lme4::lmerControl(calc.derivs = FALSE)))),
      error = function(e) NULL)
    if (is.null(fit)) next
    NC <- NC + 1L
    est <- lme4::fixef(fit); se <- sqrt(diag(as.matrix(stats::vcov(fit))))
    for (f in fnames) {
      if (!(f %in% names(est))) next
      lo <- est[f] - 1.96 * se[f]; hi <- est[f] + 1.96 * se[f]
      W[[f]] <- c(W[[f]], hi - lo)
      if (lo > rope || hi < -rope) OUT[f] <- OUT[f] + 1L
      if (lo > -rope && hi < rope) INS[f] <- INS[f] + 1L
    }
  }
  data.frame(
    param = fnames,
    true = if (!is.null(names(focal))) unname(focal) else NA_real_,
    mean_ci_width = vapply(fnames, function(f) if (length(W[[f]])) mean(W[[f]]) else NA_real_, numeric(1)),
    p_meaningful = OUT / NC, p_equivalent = INS / NC, n_converged = NC, row.names = NULL)
}

#' Sweep the number of subjects and report the ROPE decision probabilities at each size.
#' This helps identify the minimum analysable N at which a focal effect reaches a determinate
#' decision with at least some target probability (for example 0.90).
#' @export
precision_curve <- function(spec, focal, subject_ns, formula = NULL, prep = NULL, rope = 0.05, n_sims = 60) {
  if (is.character(spec)) spec <- load_spec(spec)
  if (is.null(formula)) formula <- model_formula(spec)
  if (is.null(prep)) prep <- function(d) model_data(spec, d)
  parts <- lapply(subject_ns, function(n) {
    s <- spec; s$units$subject$n <- n
    cbind(n_subject = n, precision_design(s, focal, formula = formula, prep = prep, rope = rope, n_sims = n_sims))
  })
  do.call(rbind, parts)
}
