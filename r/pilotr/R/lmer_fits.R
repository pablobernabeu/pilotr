# Shared mixed-model fitting for the replicate loops, with honest convergence accounting.
#
# The replicate loops previously wrapped each fit in suppressWarnings(), so only a hard error
# was ever visible: a boundary-singular fit and a fit whose optimiser reported non-convergence
# both counted as converged. That made `n_converged` close to meaningless for the design class
# pilotr is aimed at, because a maximal crossed random-effects structure is singular in a large
# share of replicates at realistic sample sizes (Bates et al., 2015; Matuschek et al., 2017).
#
# Singular and warning fits are still used rather than discarded. Their fixed-effect estimates
# stay interpretable, and dropping them would bias the result, since singularity is not
# independent of the variance estimates that produce it: excluding those replicates would
# preferentially remove the ones with small estimated random-effect variance, and so overstate
# precision. They are counted and reported instead, which tells the user something actionable,
# namely that the model being fitted is richer than the data can support.

# Fit one mixed model and record what the fitter actually reported. Returns the fit (NULL if it
# failed outright), whether it is boundary-singular, any warning or convergence messages, and a
# strict `converged` flag that is TRUE only when there were neither.
.fit_lmer <- function(formula, data, test = FALSE) {
  msgs <- character(0)
  # The fitter's own error message is kept rather than discarded. A model that lme4 refuses
  # outright, most often because the random-effects structure is unidentifiable at that sample
  # size, otherwise produced a result of NA with nothing to explain it, which leaves the user with
  # no way to tell an impossible model from an unlucky one.
  fit <- withCallingHandlers(
    tryCatch(
      suppressMessages(
        if (test)
          lmerTest::lmer(formula, data = data,
                         control = lme4::lmerControl(calc.derivs = FALSE))
        else
          lme4::lmer(formula, data = data,
                     control = lme4::lmerControl(calc.derivs = FALSE))),
      error = function(e) { msgs <<- c(msgs, conditionMessage(e)); NULL }),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  if (is.null(fit))
    return(list(fit = NULL, singular = FALSE, messages = msgs, converged = FALSE))
  # The optimiser records its own convergence messages separately from the R warning
  # condition, so both have to be consulted.
  opt_msgs <- tryCatch(fit@optinfo$conv$lme4$messages, error = function(e) NULL)
  singular <- isTRUE(tryCatch(lme4::isSingular(fit), error = function(e) FALSE))
  msgs <- unique(c(msgs, opt_msgs))
  list(fit = fit, singular = singular, messages = msgs,
       converged = length(msgs) == 0L && !singular)
}

# The empty per-replicate record, returned when a fit fails outright.
.fit_record_failed <- function(fnames = NULL) {
  list(fitted = FALSE, singular = FALSE, warned = FALSE, converged = FALSE)
}

# One replicate of the shared design-analysis loop: simulate, prepare, fit, and report each focal
# effect's estimate, standard error and p-value, alongside the fit diagnostics.
#
# power_mixed() and precision_design() ran separate loops that differed only in how they reduced
# the fit, and power_mixed()'s loop additionally hard-coded its own formula and data preparation
# rather than deriving them from the specification. One loop serves both, so a fix to the fitting
# or the convergence accounting reaches both at once.
#
# Kept at top level so that only the arguments travel to PSOCK workers.
#
# `test` chooses the fitter. Satterthwaite p-values come from lmerTest and cost noticeably more
# than the plain fit, so precision analysis, which needs only estimates and standard errors, asks
# for the cheaper one.
.design_rep <- function(i, spec, seeds, prep, formula, fnames, test = TRUE) {
  s <- spec; s$seed <- seeds[i]
  d <- prep(simulate_design(s, validate = FALSE))
  f <- .fit_lmer(formula, d, test = test)
  na <- stats::setNames(rep(NA_real_, length(fnames)), fnames)
  absent <- stats::setNames(logical(length(fnames)), fnames)
  if (is.null(f$fit))
    return(c(list(present = absent, est = na, se = na, p = na, coef_names = NULL,
                  error = if (length(f$messages)) f$messages[1] else NA_character_),
             .fit_record_failed()))

  co <- if (test) tryCatch(summary(f$fit)$coefficients, error = function(e) NULL) else NULL
  est <- lme4::fixef(f$fit)
  se <- sqrt(diag(as.matrix(stats::vcov(f$fit))))
  present <- stats::setNames(fnames %in% names(est), fnames)
  e <- na; s_e <- na; pv <- na
  for (fn in fnames[present]) {
    e[fn] <- est[[fn]]
    s_e[fn] <- se[[fn]]
    # lmerTest's Satterthwaite column is the p-value the power functions test against. When the
    # cheaper fitter was used, or the column is missing, the p-value stays NA and the caller
    # reports the effect as untested rather than assuming anything about it.
    if (!is.null(co) && fn %in% rownames(co) && "Pr(>|t|)" %in% colnames(co))
      pv[fn] <- co[fn, "Pr(>|t|)"]
  }
  list(present = present, est = e, se = s_e, p = pv, coef_names = names(est), fitted = TRUE,
       singular = f$singular, warned = length(f$messages) > 0L, converged = f$converged)
}

# The focal coefficient names for a specification, in the naming the auto-derived model uses.
#
# A specification's interaction keys are written "a:b", while model_data() materialises them as
# columns named "a_b", so the focal names have to follow the columns rather than the keys.
.default_focal <- function(spec) {
  vapply(names(spec$fixed$coefficients), .us, character(1), USE.NAMES = FALSE)
}

# Resolve the `focal` argument to a character vector of names and, where given, their true values.
.resolve_focal <- function(focal, spec) {
  if (is.null(focal)) {
    nms <- .default_focal(spec)
    true <- vapply(names(spec$fixed$coefficients), function(k) spec$fixed$coefficients[[k]],
                   numeric(1), USE.NAMES = FALSE)
    return(list(names = nms, true = stats::setNames(true, nms)))
  }
  if (!is.null(names(focal)) && is.numeric(focal))
    return(list(names = names(focal), true = focal))
  nms <- as.character(focal)
  list(names = nms, true = stats::setNames(rep(NA_real_, length(nms)), nms))
}

# Warn, with the fitter's own words, when not one replicate produced a fit.
#
# Every rate is then NA, which on its own gives the user nothing to act on. The usual cause is a
# random-effects structure the design cannot identify, and lme4 says so clearly, so its message is
# the most useful thing to pass on. A warning rather than an error, so that one impossible grid
# point does not abort a whole sweep.
.warn_no_fits <- function(res, n_returned, formula) {
  if (n_returned > 0L) return(invisible(NULL))
  msg <- NA_character_
  for (r in res) if (!is.null(r$error) && !is.na(r$error)) { msg <- r$error; break }
  warning(sprintf(
    "not one of the %d replicates produced a fit, so every result is NA. The model was %s.%s",
    length(res), paste(deparse(formula), collapse = " "),
    if (is.na(msg)) "" else sprintf(" The fitter reported: %s", msg)), call. = FALSE)
}

# Warn when a focal name never appeared in any fit.
#
# Such an effect yields an NA interval width and decision proportions of zero, which reads as "this
# design can decide nothing" when the real cause is a name that does not match the model. Saying so
# is the difference between a wrong answer and a question.
.warn_absent_focal <- function(seen, n_returned, coef_names) {
  missing <- names(seen)[seen == 0L & n_returned > 0L]
  if (!length(missing)) return(invisible(NULL))
  warning(sprintf(
    "focal effect%s %s never appeared in the fitted model, so %s results are empty by construction; the fitted coefficients are %s",
    if (length(missing) > 1) "s" else "", paste(sprintf("'%s'", missing), collapse = ", "),
    if (length(missing) > 1) "their" else "its",
    if (is.null(coef_names)) "unavailable" else paste(sprintf("'%s'", coef_names), collapse = ", ")),
    call. = FALSE)
}

# Collapse the per-replicate diagnostic flags into the counts reported to the user. All five are
# integers, including the two derived from arguments, so that a caller can compare them without
# tripping over the double `n_sims` came in as.
.fit_counts <- function(res, n_sims, n_returned) {
  flag <- function(nm) sum(vapply(res, function(r) isTRUE(r[[nm]]), logical(1)))
  list(n_attempted = as.integer(n_sims), n_returned = as.integer(n_returned),
       n_converged = flag("converged"), n_singular = flag("singular"),
       n_warning = flag("warned"))
}
