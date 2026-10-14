# A preregistered simplification rule, simulated replicate by replicate.
#
# power_mixed() and precision_design() fitted one model to every replicate, by default the maximal
# model the specification implies. Singular and warned fits are kept, for the reason lmer_fits.R
# gives, and that is what a researcher who reports the maximal model does. Many analysis plans
# instead name a rule: fit the maximal model, and when the fit is singular or fails to converge,
# drop the correlations, then the slopes (Barr et al., 2013; Matuschek et al., 2017). The rule
# changes both power and the Type I error rate, since the model a replicate is analysed with then
# depends on that replicate's data. The power vignette sweeps its priming design from 8 to 56
# subjects, and the maximal model was singular in at least half of the replicates at every sample
# size, so a single model could not stand in for such a plan.
#
# A ladder holds the rungs of the rule, in order, and the criterion for moving down it. Each
# replicate is fitted rung by rung and keeps the first fit the criterion accepts, or the last rung's
# fit when none is accepted, as the analysis plan would then report that model. The rung each
# replicate kept is recorded, and the result reports the share of replicates at each rung.
#
# The ladder is passed to the replicate loop as its formula. The loop may run in other processes,
# and it receives the rule together with the models the rule applies to.

# The structures model_formula() builds, which a ladder may name as its rungs.
.STRUCTURES <- c("maximal", "zcp", "intercepts")

# The ladder `formula` describes, or NULL when it is a single model: NULL, a formula or a formula
# written as text, all of which are fitted as they always were. A character vector of structure
# names gives rungs that pilotr derives with model_formula(), and the response family then chooses
# the fitter, as it does for the default model. A list gives rungs that are formulas, fitted as
# written, as a single formula given by the user is. Rungs keep their names, and take the structure
# name or their position when they have none. The rule is checked here, before any replicate runs.
.as_ladder <- function(formula, fallback, spec) {
  derived <- is.character(formula) && length(formula) > 0L &&
    !any(grepl("~", formula, fixed = TRUE))
  if (!derived && (!is.list(formula) || inherits(formula, "formula"))) return(NULL)
  if (!length(formula))
    stop("a ladder of formulas needs at least one rung", call. = FALSE)
  rung_names <- names(formula)
  if (derived) {
    unknown <- setdiff(formula, .STRUCTURES)
    if (length(unknown))
      stop(sprintf(
        "%s is not a structure model_formula() builds; name rungs among %s, or give a list of formulas",
        sprintf("'%s'", unknown[1]), paste(sprintf("'%s'", .STRUCTURES), collapse = ", ")),
        call. = FALSE)
    if (is.null(rung_names)) rung_names <- formula
    rung_names[!nzchar(rung_names)] <- formula[!nzchar(rung_names)]
    rungs <- lapply(unname(formula), function(s) model_formula(spec, s))
  } else {
    if (!all(vapply(formula, function(f) inherits(f, "formula"), logical(1))))
      stop("every rung of a list must be a formula, such as model_formula(spec, \"zcp\"); ",
           "to have pilotr derive the rungs, name their structures in a character vector, as in ",
           "c(\"maximal\", \"zcp\", \"intercepts\")", call. = FALSE)
    position <- as.character(seq_along(formula))
    if (is.null(rung_names)) rung_names <- position
    rung_names[!nzchar(rung_names)] <- position[!nzchar(rung_names)]
    rungs <- unname(formula)
  }
  twice <- unique(rung_names[duplicated(rung_names)])
  if (length(twice))
    stop(sprintf("the ladder names rung '%s' more than once; give each rung its own name",
                 twice[1]), call. = FALSE)
  structure(list(rungs = stats::setNames(rungs, rung_names), fallback = fallback,
                 derived = derived),
            class = "pilotr_ladder")
}

# Whether a rung's fit ends the descent. Under "singular", a fit that returned and is not
# boundary-singular is kept. Under "warning", it must also carry no warning or convergence message,
# which is .fit_model()'s `converged`. lme4's singular-fit notice is not among those messages, so
# the two criteria differ (see .fit_model()). A fit that failed outright is never kept while a rung
# remains.
.rung_passes <- function(f, fallback) {
  !is.null(f$fit) && !isTRUE(f$singular) && (fallback == "singular" || isTRUE(f$converged))
}

# Fit one replicate's analysis: the single model `formula`, or a ladder's rungs in turn until one
# passes the rule, else the last. A ladder's result records the rung it kept as `rung`. A rung whose
# formula is that of the rung before it, as "zcp" is for a specification whose slopes are already
# uncorrelated, reuses that fit, which is what fitting it again would return.
.fit_rule <- function(formula, data, family, test, auto) {
  if (!inherits(formula, "pilotr_ladder"))
    return(.fit_model(formula, data, family = family, test = test, auto = auto))
  rungs <- formula$rungs
  f <- NULL
  for (k in seq_along(rungs)) {
    if (k == 1L || !identical(rungs[[k]], rungs[[k - 1L]]))
      f <- .fit_model(rungs[[k]], data, family = family, test = test, auto = auto)
    if (.rung_passes(f, formula$fallback)) break
  }
  f$rung <- names(rungs)[k]
  f
}

# The formulas of an analysis, one per rung, or the single formula alone.
.rung_formulas <- function(formula) {
  if (inherits(formula, "pilotr_ladder")) formula$rungs else list(formula)
}

# The formula a replicate is left with when no fit passes: a ladder's last rung, or the single
# formula. .warn_no_fits() names it.
.last_rung <- function(formula) {
  rungs <- .rung_formulas(formula)
  rungs[[length(rungs)]]
}

# The fitter label of an analysis. Every rung of a derived ladder takes the same fitter, and a
# ladder of formulas may mix them, as when its last rung has no random terms, so each label is
# given once, in the order of the rungs.
.analysis_label <- function(formula, ...) {
  labels <- vapply(.rung_formulas(formula), .fitter_label, character(1), ...)
  paste(unique(labels), collapse = "; ")
}

# A ladder's shares as columns of a result table, one per rung, named `rung_<name>` and repeated
# on every row, since a share belongs to the run and not to one effect. The rung names are kept as
# written, where data.frame() would make them syntactic. A single model adds none.
.add_rung_columns <- function(frame, shares) {
  if (!is.null(shares)) frame[paste0("rung_", names(shares))] <- as.list(unname(shares))
  frame
}

# The share of the replicates that kept each rung, in the ladder's order. A replicate whose every
# rung failed counts at the last rung, the one the rule ended on.
.rung_shares <- function(res, ladder) {
  kept <- vapply(res, function(r) r$rung, character(1))
  counts <- tabulate(match(kept, names(ladder$rungs)), nbins = length(ladder$rungs))
  stats::setNames(counts / length(res), names(ladder$rungs))
}
