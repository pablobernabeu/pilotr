# Joint tests of several coefficients, for power_mixed().
#
# A factor with three levels enters the model as two contrast columns, and its interaction with a
# two-level factor as two product columns. The effect a study preregisters is then often the test
# of both coefficients at once, and its power differs from that of either coefficient tested by
# itself. A named list of coefficient sets in `focal` asks for that test, and each set is tested
# against the hypothesis that all its coefficients are zero.
#
# The test follows the fitter. A linear mixed model fitted by lmerTest has Satterthwaite's F from
# lmerTest::contest(), the multi-degree-of-freedom form of the t test lmerTest reports for one
# coefficient (Kuznetsova et al., 2017). A glmer() or glm() fit has the Wald chi-square of the
# set's estimates and their covariance matrix, the form of the Wald z it reports. An lm() fit has
# the F test of the same quadratic form, which equals the F test comparing the model with and
# without the set. A set of one coefficient therefore gives F = t^2, or a chi-square of z^2, and
# the p-value of the single-coefficient test. A likelihood-ratio test is not offered. It would
# refit every replicate by maximum likelihood, and for fixed effects it is anti-conservative in
# small samples (Pinheiro and Bates, 2000).

# The example every refusal below shows.
.SET_EXAMPLE <- 'list(interaction = c("A1_B", "A2_B"))'

# A list `focal` checked as a set of named, non-empty sets of distinct coefficient names. A set is
# matched to its results by its name, so an unnamed or duplicated one is refused.
.as_focal_sets <- function(focal) {
  nms <- names(focal)
  if (!length(focal) || is.null(nms) || anyNA(nms) || any(!nzchar(nms)) || anyDuplicated(nms))
    stop("a list `focal` holds sets of coefficients, each tested jointly, and every set needs ",
         "its own name, as in ", .SET_EXAMPLE, call. = FALSE)
  for (k in nms) {
    s <- focal[[k]]
    if (!is.character(s) || !length(s) || anyNA(s) || any(!nzchar(s)) || anyDuplicated(s))
      stop(sprintf("set '%s' of `focal` must name one or more distinct coefficients, as in %s",
                   k, .SET_EXAMPLE), call. = FALSE)
  }
  lapply(focal, as.character)
}

# Checks made before any worker starts: the sets are well formed, and the test is two-sided. A
# joint test asks whether any coefficient of the set differs from zero, in whatever direction.
.check_focal_sets <- function(focal, alternative) {
  if (!is.list(focal)) return(invisible(NULL))
  .as_focal_sets(focal)
  if (alternative != "two.sided")
    stop("a set of coefficients is tested jointly, by an F or chi-square test that has no ",
         "direction, so `alternative` must be \"two.sided\" when `focal` is a list. For a ",
         "one-sided test, name a single coefficient.", call. = FALSE)
  invisible(NULL)
}

# precision_design() and precision_curve() build an interval around each coefficient's estimate,
# and a set of coefficients has no single estimate to build one around.
.refuse_focal_sets <- function(focal, fn) {
  if (is.list(focal))
    stop(sprintf(paste0(
      "%s builds an interval around each coefficient's estimate, and a set of coefficients has ",
      "no single estimate. Name the coefficients one by one, or test the set jointly with ",
      "power_mixed()."), fn), call. = FALSE)
  invisible(NULL)
}

# The joint test that the coefficients `set` of `fit` are all zero, as a list of its statistic, its
# degrees of freedom and its p-value. `df2` is NA for a chi-square. The result is NULL when a
# coefficient of the set is missing from the fit or was not estimated, or when the test fails. It
# is also NULL for a plain lme4::lmer() fit, which lacks lmerTest's degrees of freedom and leaves
# a single coefficient untested too.
.joint_test <- function(fit, set) {
  b <- if (inherits(fit, "merMod")) lme4::fixef(fit) else stats::coef(fit)
  if (!all(set %in% names(b)) || anyNA(b[set])) return(NULL)
  tryCatch({
    if (inherits(fit, "lmerModLmerTest")) {
      # One row of L per coefficient of the set, picking it out of the fixed effects.
      L <- matrix(0, length(set), length(b))
      L[cbind(seq_along(set), match(set, names(b)))] <- 1
      ct <- lmerTest::contest(fit, L, joint = TRUE)
      return(list(statistic = ct[["F value"]], df1 = ct[["NumDF"]], df2 = ct[["DenDF"]],
                  p = ct[["Pr(>F)"]]))
    }
    q <- as.numeric(length(set))
    V <- as.matrix(stats::vcov(fit))[set, set, drop = FALSE]
    w <- sum(b[set] * solve(V, b[set]))
    if (inherits(fit, "glmerMod") || inherits(fit, "glm"))
      list(statistic = w, df1 = q, df2 = NA_real_, p = stats::pchisq(w, q, lower.tail = FALSE))
    else if (inherits(fit, "lm")) {
      df2 <- stats::df.residual(fit)
      list(statistic = w / q, df1 = q, df2 = df2,
           p = stats::pf(w / q, q, df2, lower.tail = FALSE))
    }
  }, error = function(e) NULL)
}

# The record .design_rep() returns when `focal` holds sets: each set's joint p-value, with no
# estimate or standard error, since a set has no single one, and the fit's diagnostics, as for a
# single coefficient.
.joint_record <- function(f, sets) {
  coefs <- if (inherits(f$fit, "merMod")) names(lme4::fixef(f$fit)) else names(stats::coef(f$fit))
  na <- stats::setNames(rep(NA_real_, length(sets)), names(sets))
  p <- vapply(sets, function(s) {
    t <- .joint_test(f$fit, s)
    if (is.null(t)) NA_real_ else t$p
  }, numeric(1))
  list(present = vapply(sets, function(s) all(s %in% coefs), logical(1)),
       est = na, se = na, p = p, coef_names = coefs, fitted = TRUE,
       singular = f$singular, warned = length(f$messages) > 0L, converged = f$converged)
}

# Warn when a set was never tested, naming the coefficients the fitted model lacks. Its results are
# then NA, and the warning says why, as .warn_absent_focal() does for a single coefficient.
.warn_untested_sets <- function(sets, seen, n_returned, coef_names) {
  never <- names(seen)[seen == 0L & n_returned > 0L]
  if (!length(never)) return(invisible(NULL))
  quote <- function(x) paste(sprintf("'%s'", x), collapse = ", ")
  lacking <- unique(unlist(lapply(sets[never], setdiff, coef_names), use.names = FALSE))
  several <- length(never) > 1L
  warning(sprintf(
    "focal set%s %s %s never tested, since %s, so %s results are empty by construction. The fitted coefficients are %s.",
    if (several) "s" else "", quote(never), if (several) "were" else "was",
    if (length(lacking)) sprintf("the fitted model lacks %s", quote(lacking))
    else "the joint test failed in every replicate",
    if (several) "their" else "its",
    if (is.null(coef_names)) "unavailable" else quote(coef_names)),
    call. = FALSE)
}
