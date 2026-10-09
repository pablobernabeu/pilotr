# Shared mixed-model fitting for the replicate loops, with honest convergence accounting.
#
# The replicate loops previously wrapped each fit in suppressWarnings(), so only a hard error
# was ever visible: a boundary-singular fit and a fit whose optimiser reported non-convergence
# both counted as converged. That made `n_converged` close to meaningless for the design class
# pilotr is aimed at, because a maximal crossed random-effects structure is singular in a large
# share of replicates at realistic sample sizes (Bates et al., 2015; Matuschek et al., 2017).
#
# Singular and warning fits still count towards the result. Their fixed-effect estimates
# stay interpretable, and dropping them would bias the result, since singularity is not
# independent of the variance estimates that produce it: excluding those replicates would
# preferentially remove the ones with small estimated random-effect variance, and so overstate
# precision. They are counted and reported instead, which tells the user something actionable,
# namely that the model being fitted is richer than the data can support.

# The fitter follows the response family. Earlier versions fitted every replicate by lmer(),
# whatever the family. A bernoulli or poisson specification writes its coefficients on the logit
# or log scale, and a linear model of the 0/1 or count response estimates a difference on the
# response scale instead. Compared with the true value, that difference made Type M read 0.21 on a
# crossed accuracy design that a logistic model shows exaggerating, at 1.12, and set ROPE intervals
# on one scale against a region on another. Power hardly moved, since the two tests agree on the
# sign and nearly always on significance. glmer() fits those two families on their own link scale,
# with a Wald z test. A model with no random terms, such as that of every shipped
# between-subjects example, made lmer() stop with "No random effects terms specified in formula"
# in every replicate. Such a model is fitted by lm() or glm().
#
# Ordinal and beta responses have no frequentist model here yet, since a proportional-odds or
# Beta mixed model would need a package beyond lme4. They are fitted by the linear fallback on
# the response scale. The callers withhold the quantities that compare that scale with the
# specification's link scale, and raise the warnings below.

# The families whose coefficients are on a link scale, with the scale's name.
.LINK_SCALE <- c(bernoulli = "logit", poisson = "log", ordinal = "logit", beta = "logit")

# The two of those that lme4 and stats fit on that scale, with the family object's name.
.GLM_FAMILY <- c(bernoulli = "binomial", poisson = "poisson")

# The warnings for a family fitted by the linear fallback. power_mixed() in the Python twin raises
# the first word for word, for every family in .LINK_SCALE, since it has no other model.
.LINEAR_POWER <- paste0(
  "power_mixed() fits a linear model to the %s response on its own scale, while the ",
  "specification's coefficients are on the %s scale, so the mean estimate and Type M are ",
  "withheld (NA) and Type S compares signs only. For a model on the link scale, use ",
  "generate_design_analysis() or brms_bridge() in the R package.")
.LINEAR_PRECISION <- paste0(
  "precision_design() fits a linear model to the %s response on its own scale, while the ",
  "region of practical equivalence and the specification's coefficients are on the %s scale, ",
  "so the decision probabilities are withheld (NA). For a decision on the link scale, use ",
  "generate_design_analysis() in the R package.")

# Whether a model formula has a random-effects term, a `|` or `||` call anywhere on its right-hand
# side, which is where lme4 looks for one. lme4::findbars() answers the same question, but lme4
# 2.0 moved it to the reformulas package and warns when it is called through lme4. The fitter
# needs to know only whether there is such a term. lmer() also takes a formula written as text,
# which is read as a formula first. Taken as it came, the text showed no such term, and lm() then
# fitted the model without its random effects.
.has_bars <- function(formula) {
  formula <- stats::as.formula(formula)
  walk <- function(x) is.call(x) &&
    (identical(x[[1L]], as.name("|")) || identical(x[[1L]], as.name("||")) ||
       any(vapply(as.list(x)[-1L], walk, logical(1))))
  walk(formula[[length(formula)]])
}

# Which model a fit takes: the GLM family name (NULL for a linear model) and whether the formula
# has random terms. The family chooses only for the model pilotr derives (`auto`). A formula the
# user writes is fitted as a linear model, as it always was, since only its author knows the
# scale its response is on.
.model_kind <- function(formula, family, auto) {
  glm_family <- if (isTRUE(auto) && family %in% names(.GLM_FAMILY)) .GLM_FAMILY[[family]]
  list(glm_family = glm_family, has_re = .has_bars(formula))
}

# Whether a call fits the linear fallback to a response whose coefficients are on a link scale,
# which is when the scale-dependent quantities are withheld.
.linear_on_link <- function(family, auto) {
  isTRUE(auto) && family %in% setdiff(names(.LINK_SCALE), names(.GLM_FAMILY))
}

# The fitter and the test behind a run's estimates, as the result reports them, for example
# "lme4::glmer (binomial), Wald z". With `test = FALSE`, the only inference made is the interval
# precision_design() builds, the estimate plus or minus 1.96 standard errors, a Wald z interval.
.fitter_label <- function(formula, family, test, auto) {
  k <- .model_kind(formula, family, auto)
  fitter <- if (!is.null(k$glm_family))
    sprintf("%s (%s)", if (k$has_re) "lme4::glmer" else "stats::glm", k$glm_family)
  else if (k$has_re) { if (test) "lmerTest::lmer" else "lme4::lmer" }
  else "stats::lm"
  inference <- if (!test || !is.null(k$glm_family)) "Wald z"
    else if (k$has_re) "Satterthwaite t" else "t"
  paste(fitter, inference, sep = ", ")
}

# Fit one replicate's model and record what the fitter reported. Returns the fit (NULL if it
# failed outright), whether it is boundary-singular, any warning or convergence messages other than
# the singular-fit notice, and a strict `converged` flag that is TRUE only when there were neither.
# `family` is the specification's response family, and `auto` says whether pilotr derived the
# formula. Together, they choose the fitter, as .model_kind() describes.
.fit_model <- function(formula, data, family = "gaussian", test = FALSE, auto = TRUE) {
  k <- .model_kind(formula, family, auto)
  fam <- if (!is.null(k$glm_family)) switch(k$glm_family, binomial = stats::binomial(),
                                            poisson = stats::poisson())
  msgs <- character(0)
  # By default lme4 stores its boundary (singular) fit notice with the optimiser's messages, which
  # are read below. Every singular fit was then also counted as a fit with warnings, and n_warning
  # could never differ from n_singular. The check is switched off here, and singularity is taken
  # from isSingular(), whose default tolerance matches this one. Filtering the notice by its text
  # would break whenever lme4 rewords it, as it did in 1.1-21. glmer() takes the same settings.
  singular_cc <- lme4::.makeCC(action = "ignore", tol = 1e-4)
  # The fitter's own error message is kept and passed on. A model that lme4 refuses
  # outright, most often because the random-effects structure is unidentifiable at that sample
  # size, otherwise produced a result of NA with nothing to explain it, which leaves the user with
  # no way to tell an impossible model from an unlucky one.
  fit <- withCallingHandlers(
    tryCatch(
      suppressMessages(
        if (k$has_re && !is.null(fam))
          lme4::glmer(formula, data = data, family = fam,
                      control = lme4::glmerControl(calc.derivs = FALSE,
                                                   check.conv.singular = singular_cc))
        else if (k$has_re && test)
          lmerTest::lmer(formula, data = data,
                         control = lme4::lmerControl(calc.derivs = FALSE,
                                                     check.conv.singular = singular_cc))
        else if (k$has_re)
          lme4::lmer(formula, data = data,
                     control = lme4::lmerControl(calc.derivs = FALSE,
                                                 check.conv.singular = singular_cc))
        else if (!is.null(fam))
          stats::glm(formula, data = data, family = fam)
        else
          stats::lm(formula, data = data)),
      error = function(e) { msgs <<- c(msgs, conditionMessage(e)); NULL }),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  if (is.null(fit))
    return(list(fit = NULL, singular = FALSE, messages = msgs, converged = FALSE))
  if (inherits(fit, "merMod")) {
    # The optimiser records its own convergence messages separately from the R warning
    # condition, so both have to be consulted.
    opt_msgs <- tryCatch(fit@optinfo$conv$lme4$messages, error = function(e) NULL)
    singular <- isTRUE(tryCatch(lme4::isSingular(fit), error = function(e) FALSE))
  } else {
    # lm() and glm() estimate no variance component, so neither fit can be singular. glm() also
    # records in the fit whether its iterations converged, beside the warning it raises. The
    # message below is glm.fit's own, so that unique() keeps one of the two.
    opt_msgs <- if (inherits(fit, "glm") && !isTRUE(fit$converged))
      "glm.fit: algorithm did not converge"
    singular <- FALSE
  }
  msgs <- unique(c(msgs, opt_msgs))
  list(fit = fit, singular = singular, messages = msgs,
       converged = length(msgs) == 0L && !singular)
}

# A model formula fitted by lmer(), or by lm() when it has no random terms, whatever the family.
.fit_lmer <- function(formula, data, test = FALSE) {
  .fit_model(formula, data, family = "gaussian", test = test, auto = FALSE)
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
# instead of deriving them from the specification. One loop serves both, so a fix to the fitting
# or the convergence accounting reaches both at once.
#
# Kept at top level so that only the arguments travel to PSOCK workers.
#
# `test` chooses between the linear fitters. Satterthwaite p-values come from lmerTest and cost
# noticeably more than the plain fit, so precision analysis, which needs only estimates and standard
# errors, asks for the cheaper one. `auto` says whether pilotr derived the formula, which lets the
# response family choose the fitter (see .fit_model()).
.design_rep <- function(i, spec, seeds, prep, formula, fnames, test = TRUE, auto = TRUE) {
  s <- spec; s[["seed"]] <- seeds[i]
  d <- prep(simulate_design(s, validate = FALSE))
  f <- .fit_model(formula, d, family = spec[["response"]][["family"]], test = test, auto = auto)
  na <- stats::setNames(rep(NA_real_, length(fnames)), fnames)
  absent <- stats::setNames(logical(length(fnames)), fnames)
  if (is.null(f$fit))
    return(c(list(present = absent, est = na, se = na, p = na, coef_names = NULL,
                  error = if (length(f$messages)) f$messages[1] else NA_character_),
             .fit_record_failed()))

  co <- if (test) tryCatch(summary(f$fit)$coefficients, error = function(e) NULL) else NULL
  # The test's column is named for its reference distribution: Pr(>|t|) for lmerTest's
  # Satterthwaite test and for lm(), Pr(>|z|) for the Wald z of glm() and glmer().
  pcol <- if (is.null(co)) NA_character_ else intersect(c("Pr(>|t|)", "Pr(>|z|)"), colnames(co))[1]
  est <- if (inherits(f$fit, "merMod")) lme4::fixef(f$fit) else stats::coef(f$fit)
  se <- sqrt(diag(as.matrix(stats::vcov(f$fit))))
  present <- stats::setNames(fnames %in% names(est), fnames)
  e <- na; s_e <- na; pv <- na
  for (fn in fnames[present]) {
    e[fn] <- est[[fn]]
    s_e[fn] <- se[[fn]]
    # That column holds the p-value the power functions test against. When the cheaper fitter was
    # used, or the column is missing, the p-value stays NA and the caller reports the effect as
    # untested, and assumes nothing about it. A coefficient that lm() could not estimate, whose
    # row summary() leaves out, is treated in the same way.
    if (!is.na(pcol) && fn %in% rownames(co))
      pv[fn] <- co[fn, pcol]
  }
  list(present = present, est = e, se = s_e, p = pv, coef_names = names(est), fitted = TRUE,
       singular = f$singular, warned = length(f$messages) > 0L, converged = f$converged)
}

# The focal coefficient names for a specification, in the naming the auto-derived model uses.
#
# A specification's interaction keys are written "a:b", while model_data() materialises them as
# columns named "a_b", so the focal names have to follow the columns.
.default_focal <- function(spec) {
  vapply(names(spec[["fixed"]][["coefficients"]]), .us, character(1), USE.NAMES = FALSE)
}

# Resolve the `focal` argument to a character vector of names and, where given, their true values.
.resolve_focal <- function(focal, spec) {
  if (is.null(focal)) {
    nms <- .default_focal(spec)
    coeffs <- spec[["fixed"]][["coefficients"]]
    true <- vapply(names(coeffs), function(k) coeffs[[k]], numeric(1), USE.NAMES = FALSE)
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
# the most useful thing to pass on. It is raised as a warning, so that one impossible grid
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
