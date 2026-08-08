# Specification validation and version negotiation.
#
# load_spec() was a bare jsonlite::fromJSON call: no validation and no version check. A strict
# draft-07 schema has shipped at spec/design.schema.json since 0.1, and no code path consulted
# it, so a specification with a misspelled field loaded silently and kept the misspelling.
#
# That matters more than a typo usually would, because several of the ways a specification can
# be wrong produce plausible data rather than an error:
#
#   * A mistyped coefficient key resolves to no column, so that effect is silently set to zero.
#     A spec whose focal effect is named "cnod" instead of "cond" generates exactly the data of
#     a null design, and reports success.
#   * varies_by took anything other than "subject" as item-level, so a per-trial predictor was
#     silently given one value per item.
#   * A response parameter belonging to another family is ignored rather than refused.
#
# It also matters for version negotiation. Once 0.3 features exist, a 0.3 specification opened
# by a 0.2 implementation, an un-upgraded Python twin, or a cached browser build produces
# different and wrong data while reporting success. A 0.2 implementation has no version check
# and cannot be fixed retrospectively, but a specification that uses a 0.3 feature can be made
# to say so, and every implementation from 0.3 onwards refuses what it does not understand.

# The specification version this implementation writes and understands.
.SPEC_VERSION <- "0.3"

# Features introduced in 0.3. A specification using any of them is read differently by a 0.2
# implementation, so it has to declare 0.3 or later.
.spec_0_3_features <- function(spec) {
  found <- character(0)
  for (p in .orelse(spec$predictors, list())) {
    if (identical(p$varies_by, "observation")) found <- c(found, 'predictors varies_by "observation"')
    if (!is.null(p$dist)) found <- c(found, "predictors dist")
    if (!is.null(p$reliability) && !identical(as.numeric(p$reliability), 1))
      found <- c(found, "predictors reliability")
  }
  if (identical(spec$response$family, "exgaussian")) found <- c(found, 'the "exgaussian" family')
  for (g in names(.orelse(spec$random, list()))) {
    if (!is.null(spec$random[[g]]$correlated)) found <- c(found, "random correlated")
    if (any(grepl(":", names(spec$random[[g]]$slopes), fixed = TRUE)))
      found <- c(found, "interaction random slopes")
  }
  unique(found)
}

# Not `%||%`: base R gained that operator in 4.4.0, above the 4.0.0 pilotr declares, so defining
# it here would shadow the base version on new R and be the only definition on old R.
.orelse <- function(a, b) if (is.null(a)) b else a

# Resolve a specification argument to a validated list. Every public entry point calls this
# exactly once, so that a path is read once and validation runs once for a whole replicate loop
# rather than once per replicate.
.as_spec <- function(spec, strict = TRUE) {
  if (is.character(spec)) return(load_spec(spec, validate = strict))
  validate_spec(spec, strict = strict)
  spec
}

.is_scalar_string <- function(x) is.character(x) && length(x) == 1L && !is.na(x)
.is_scalar_number <- function(x) is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x)
.is_whole <- function(x) .is_scalar_number(x) && x == round(x)

# A version arrives as whatever JSON produced. A single part is read as a whole version, so "1"
# and the JSON number 1.0 both mean 1.0. The padding is what keeps the two engines agreeing about
# the same file: R renders the number 1.0 as "1" and Python as "1.0", so without it one
# implementation called the specification malformed while the other read it as version 1.
.parse_version <- function(v) {
  s <- if (is.numeric(v)) format(v, digits = 15) else as.character(v)
  parts <- suppressWarnings(as.integer(strsplit(s, ".", fixed = TRUE)[[1]]))
  if (length(parts) == 1L) parts <- c(parts, 0L)
  if (length(parts) < 2L || anyNA(parts[1:2])) return(NULL)
  parts[1:2]
}

# Response families and the parameters each one uses. Anything else supplied under `response`
# is refused rather than ignored, because a leftover parameter from another family is usually a
# half-finished edit and silently dropping it hides the mistake.
.family_params <- list(
  gaussian          = "sigma",
  lognormal         = "sigma",
  shifted_lognormal = c("sigma", "shift"),
  exgaussian        = c("sigma", "beta"),
  bernoulli         = character(0),
  poisson           = character(0),
  ordinal           = "thresholds",
  beta              = "phi"
)

# Families whose response value is rounded when `response.round` is set. For the others the
# outcome is an integer already, so `round` would do nothing and is refused as a likely mistake.
.rounding_families <- c("gaussian", "lognormal", "shifted_lognormal", "exgaussian", "beta")

#' Validate a design specification
#'
#' Check a design specification against the portable schema and against the cross-field rules
#' the schema cannot express, and check that its declared `spec_version` is one this
#' implementation understands. Called by [load_spec()] by default.
#'
#' @details
#' Validation exists because several ways of getting a specification wrong produce plausible
#' data rather than an error. A mistyped coefficient key resolves to no column and so silently
#' sets that effect to zero, which generates exactly the data of a null design and reports
#' success. A response parameter left over from another family is ignored. Neither is
#' detectable in the output, which is why they are refused here.
#'
#' Version negotiation covers the other direction. A specification that uses a feature
#' introduced in 0.3 is read differently by a 0.2 implementation, so it must declare 0.3 or
#' later; a specification declaring a version newer than this implementation is refused rather
#' than partially understood. A specification with no `spec_version` is treated as 0.2, which is
#' what every specification written before the field existed is.
#'
#' @param spec A design specification (path or list).
#' @param strict Whether an unrecognised field is an error (the default) or a warning. Set
#'   `FALSE` to load a specification carrying private annotations, accepting that a misspelled
#'   field will then be ignored rather than reported.
#' @return The specification, invisibly, so that the call can be chained.
#' @examples
#' spec <- build_spec(list(name = "demo", seed = 1, design_kind = "between",
#'   factor_name = "group", lev1 = "a", lev2 = "b", n_subject = 20,
#'   intercept = 0, effect = 0.5, family = "gaussian", resp_name = "", sigma = 1))
#' validate_spec(spec)
#'
#' # A mistyped coefficient key is refused rather than silently treated as a zero effect.
#' bad <- spec
#' bad$fixed$coefficients <- list(effct = 0.5)
#' try(validate_spec(bad))
#' @export
validate_spec <- function(spec, strict = TRUE) {
  if (is.character(spec)) spec <- load_spec(spec, validate = FALSE)
  if (!is.list(spec) || is.null(names(spec)))
    stop("a design specification must be a named list (a JSON object)", call. = FALSE)

  problems <- character(0)
  soft <- character(0)
  bad <- function(...) problems <<- c(problems, paste0(...))
  unknown <- function(...) if (strict) bad(...) else soft <<- c(soft, paste0(...))

  # ---- version ----
  declared <- .orelse(spec$spec_version, "0.2")
  dv <- .parse_version(declared)
  sv <- .parse_version(.SPEC_VERSION)
  if (is.null(dv)) {
    bad("spec_version '", declared, "' is not of the form 'major.minor'")
  } else {
    # The version as pilotr read it, rather than as it was written, so that the two engines report
    # the same thing about a JSON number they render differently.
    shown <- paste0(dv[1], ".", dv[2])
    if (dv[1] > sv[1] || (dv[1] == sv[1] && dv[2] > sv[2]))
      bad("this specification declares spec_version ", shown,
          ", which is newer than the ", .SPEC_VERSION,
          " this version of pilotr understands; please upgrade pilotr")
    used <- .spec_0_3_features(spec)
    if (length(used) && (dv[1] == 0L && dv[2] < 3L))
      bad("this specification uses ", paste(used, collapse = ", "),
          ", which requires spec_version \"0.3\", but declares ", shown,
          "; a 0.2 implementation would read it differently and silently generate different data")
  }

  # ---- top level ----
  known_top <- c("spec_version", "name", "seed", "units", "factors", "predictors",
                 "fixed", "random", "response")
  for (k in setdiff(names(spec), known_top))
    unknown("unknown top-level field '", k, "'; expected one of ", paste(known_top, collapse = ", "))
  for (k in c("name", "seed", "units", "fixed", "response"))
    if (is.null(spec[[k]])) bad("required top-level field '", k, "' is missing")

  if (!is.null(spec$name) && !.is_scalar_string(spec$name)) bad("'name' must be a single string")
  if (!is.null(spec$seed) && !.is_whole(spec$seed)) bad("'seed' must be a single whole number")

  # ---- units ----
  u <- spec$units
  has_item <- FALSE
  if (!is.null(u)) {
    if (!is.list(u) || is.null(names(u))) bad("'units' must be an object")
    else {
      for (k in setdiff(names(u), c("subject", "item")))
        unknown("unknown unit '", k, "'; only 'subject' and 'item' exist")
      if (is.null(u$subject)) bad("'units.subject' is required")
      has_item <- !is.null(u$item)
      for (nm in intersect(names(u), c("subject", "item"))) {
        un <- u[[nm]]
        # Only the list test, not the names() test used elsewhere for an object: an empty JSON
        # object arrives as an unnamed empty list, and the twin reports the missing n for it
        # rather than the wrong shape. Without this guard `un$n` below threw the base error
        # "$ operator is invalid for atomic vectors" where Python said what was wrong.
        if (!is.list(un)) { bad("'units.", nm, "' must be an object"); next }
        for (k in setdiff(names(un), c("n", "per_subject")))
          unknown("unknown field 'units.", nm, ".", k, "'")
        if (!.is_whole(un$n) || un$n < 1) bad("'units.", nm, ".n' must be a whole number of at least 1")
        if (!is.null(un$per_subject)) {
          if (!identical(nm, "item"))
            bad("'per_subject' belongs to 'units.item', not 'units.", nm, "'")
          else if (!.is_whole(un$per_subject) || un$per_subject < 1)
            bad("'units.item.per_subject' must be a whole number of at least 1")
          else if (.is_whole(un$n) && un$per_subject > un$n)
            bad("'units.item.per_subject' (", un$per_subject,
                ") cannot exceed the number of items (", un$n, ")")
        }
      }
    }
  }

  # ---- factors ----
  contrast_cols <- character(0)
  if (!is.null(spec$factors)) {
    if (!is.list(spec$factors) || !is.null(names(spec$factors)))
      bad("'factors' must be an array of factor objects")
    else for (i in seq_along(spec$factors)) {
      f <- spec$factors[[i]]; where <- paste0("factors[", i, "]")
      if (!is.list(f) || is.null(names(f))) { bad(where, " must be an object"); next }
      for (k in setdiff(names(f), c("name", "levels", "contrasts", "vary_within", "between")))
        unknown("unknown field '", where, ".", k, "'")
      if (!.is_scalar_string(f$name)) bad(where, ".name must be a single string")
      nlev <- length(f$levels)
      if (!is.character(f$levels) || nlev < 2)
        bad(where, ".levels must be an array of at least two strings")
      if (!is.list(f$contrasts) || is.null(names(f$contrasts)) || !length(f$contrasts))
        bad(where, ".contrasts must be a non-empty object mapping contrast columns to one value per level")
      else for (cn in names(f$contrasts)) {
        contrast_cols <- c(contrast_cols, cn)
        v <- f$contrasts[[cn]]
        if (!is.numeric(v) || anyNA(v)) bad(where, ".contrasts.", cn, " must be numeric")
        else if (nlev >= 2 && length(v) != nlev)
          bad(where, ".contrasts.", cn, " has ", length(v), " value(s) but the factor has ",
              nlev, " level(s)")
      }
      # A single string is accepted where an array belongs. pilotr's own spec_json() emitted that
      # form before 0.3, because a blanket auto_unbox collapsed every one-element array, so
      # refusing it would mean refusing files pilotr itself wrote. The reading is unambiguous and
      # both engines already treat the two alike.
      if (!is.null(f$vary_within)) {
        if (!is.character(f$vary_within) || !length(f$vary_within))
          bad(where, ".vary_within must be a unit name or an array of unit names")
        else for (w in f$vary_within) {
          if (!w %in% c("subject", "item")) bad(where, ".vary_within contains '", w,
                                               "'; only 'subject' and 'item' are allowed")
          else if (identical(w, "item") && !has_item)
            bad(where, ".vary_within names 'item' but the design has no item unit")
        }
      }
      if (!is.null(f$between)) {
        if (!.is_scalar_string(f$between) || !f$between %in% c("subject", "item"))
          bad(where, ".between must be 'subject' or 'item'")
        else if (identical(f$between, "item") && !has_item)
          bad(where, ".between is 'item' but the design has no item unit")
      }
      if (is.null(f$vary_within) && is.null(f$between))
        bad(where, " must set either 'vary_within' or 'between'")
    }
  }

  # ---- predictors ----
  pred_names <- character(0)
  if (!is.null(spec$predictors)) {
    if (!is.list(spec$predictors) || !is.null(names(spec$predictors)))
      bad("'predictors' must be an array of predictor objects")
    else for (i in seq_along(spec$predictors)) {
      p <- spec$predictors[[i]]; where <- paste0("predictors[", i, "]")
      if (!is.list(p) || is.null(names(p))) { bad(where, " must be an object"); next }
      for (k in setdiff(names(p), c("name", "varies_by", "mean", "sd", "dist",
                                    "min", "max", "reliability")))
        unknown("unknown field '", where, ".", k, "'")
      if (!.is_scalar_string(p$name)) bad(where, ".name must be a single string")
      else pred_names <- c(pred_names, p$name)
      if (!.is_scalar_string(p$varies_by) ||
          !p$varies_by %in% c("subject", "item", "observation"))
        bad(where, ".varies_by must be 'subject', 'item' or 'observation'",
            if (!is.null(p$varies_by)) paste0(", not '", p$varies_by, "'") else "")
      else if (identical(p$varies_by, "item") && !has_item)
        bad(where, ".varies_by is 'item' but the design has no item unit")
      dist <- .orelse(p$dist, "normal")
      if (!.is_scalar_string(dist) || !dist %in% c("normal", "uniform"))
        bad(where, ".dist must be 'normal' or 'uniform'")
      else if (identical(dist, "uniform")) {
        if (!.is_scalar_number(p$min) || !.is_scalar_number(p$max))
          bad(where, " uses dist 'uniform' and so needs numeric 'min' and 'max'")
        else if (p$min >= p$max) bad(where, ".min must be less than ", where, ".max")
        for (k in intersect(names(p), c("mean", "sd")))
          unknown(where, ".", k, " is ignored when dist is 'uniform'")
      } else {
        for (k in intersect(names(p), c("min", "max")))
          unknown(where, ".", k, " is ignored when dist is 'normal'")
        if (!is.null(p$mean) && !.is_scalar_number(p$mean)) bad(where, ".mean must be a number")
        if (!is.null(p$sd) && (!.is_scalar_number(p$sd) || p$sd < 0))
          bad(where, ".sd must be a number of at least 0")
      }
      if (!is.null(p$reliability) &&
          (!.is_scalar_number(p$reliability) || p$reliability <= 0 || p$reliability > 1))
        bad(where, ".reliability must be greater than 0 and at most 1")
    }
  }
  if (anyDuplicated(pred_names)) bad("duplicated predictor name(s): ",
                                     paste(unique(pred_names[duplicated(pred_names)]), collapse = ", "))
  known_cols <- c(contrast_cols, pred_names)

  # Every coefficient and slope key must resolve to a contrast column or a predictor. An
  # unresolved key contributes zero, so a typo silently removes the effect.
  check_key <- function(key, where) {
    parts <- strsplit(key, ":", fixed = TRUE)[[1]]
    miss <- setdiff(parts, known_cols)
    if (length(miss))
      bad(where, " '", key, "' names ", paste(sprintf("'%s'", miss), collapse = ", "),
          ", which ", if (length(miss) > 1) "are" else "is",
          " neither a contrast column nor a predictor; available columns are ",
          if (length(known_cols)) paste(sprintf("'%s'", known_cols), collapse = ", ") else "(none)",
          ". An unresolved key contributes zero, so this would silently drop the term")
  }

  # ---- fixed ----
  fx <- spec$fixed
  if (!is.null(fx)) {
    if (!is.list(fx) || is.null(names(fx))) bad("'fixed' must be an object")
    else {
      for (k in setdiff(names(fx), c("intercept", "coefficients")))
        unknown("unknown field 'fixed.", k, "'")
      if (!.is_scalar_number(fx$intercept)) bad("'fixed.intercept' must be a single number")
      if (is.null(fx$coefficients)) bad("'fixed.coefficients' is required (use {} for none)")
      else if (!is.list(fx$coefficients)) bad("'fixed.coefficients' must be an object")
      else for (k in names(fx$coefficients)) {
        if (!.is_scalar_number(fx$coefficients[[k]]))
          bad("'fixed.coefficients.", k, "' must be a single number")
        check_key(k, "fixed.coefficients")
      }
    }
  }

  # ---- random ----
  if (!is.null(spec$random) && length(spec$random)) {
    if (!is.list(spec$random) || is.null(names(spec$random)))
      bad("'random' must be an object keyed by grouping factor")
    else for (g in names(spec$random)) {
      re <- spec$random[[g]]; where <- paste0("random.", g)
      if (!is.list(re) || is.null(names(re))) { bad(where, " must be an object"); next }
      for (k in setdiff(names(re), c("intercept_sd", "slopes", "correlations", "correlated",
                                     "over", "n")))
        unknown("unknown field '", where, ".", k, "'")
      if (!.is_scalar_number(re$intercept_sd) || re$intercept_sd < 0)
        bad(where, ".intercept_sd is required and must be at least 0")
      cols <- c("intercept", names(re$slopes))
      if (!is.null(re$slopes)) {
        if (!is.list(re$slopes)) bad(where, ".slopes must be an object")
        else for (k in names(re$slopes)) {
          if (!.is_scalar_number(re$slopes[[k]]) || re$slopes[[k]] < 0)
            bad(where, ".slopes.", k, " must be a number of at least 0")
          check_key(k, paste0(where, ".slopes"))
        }
      }
      if (!is.null(re$correlations)) {
        if (!is.list(re$correlations)) bad(where, ".correlations must be an object")
        else for (k in names(re$correlations)) {
          v <- re$correlations[[k]]
          if (!.is_scalar_number(v) || v < -1 || v > 1)
            bad(where, ".correlations.", k, " must be between -1 and 1")
          parts <- trimws(strsplit(gsub("~", ",", k), ",")[[1]])
          if (length(parts) != 2L)
            bad(where, ".correlations key '", k, "' must name two terms, as 'a,b'")
          else {
            miss <- setdiff(parts, cols)
            if (length(miss))
              bad(where, ".correlations key '", k, "' names ",
                  paste(sprintf("'%s'", miss), collapse = ", "),
                  ", which is not a random-effect term of ", g,
                  "; its terms are ", paste(sprintf("'%s'", cols), collapse = ", "))
          }
        }
      }
      if (!is.null(re$correlated) && !(is.logical(re$correlated) && length(re$correlated) == 1L))
        bad(where, ".correlated must be TRUE or FALSE")
      if (isTRUE(identical(re$correlated, FALSE)) && length(re$correlations))
        bad(where, " sets correlated = false but also supplies correlations; one of the two has to go")
      if (g %in% c("subject", "item")) {
        for (k in intersect(names(re), c("over", "n")))
          bad(where, ".", k, " applies only to an extra grouping factor, not to '", g, "'")
      } else {
        if (!.is_scalar_string(re$over) || !re$over %in% c("subject", "item"))
          bad(where, ".over is required for an extra grouping factor and must be 'subject' or 'item'")
        else if (identical(re$over, "item") && !has_item)
          bad(where, ".over is 'item' but the design has no item unit")
        if (!.is_whole(re$n) || re$n < 1)
          bad(where, ".n is required for an extra grouping factor and must be a whole number of at least 1")
      }
    }
  }

  # ---- response ----
  r <- spec$response
  if (!is.null(r)) {
    if (!is.list(r) || is.null(names(r))) bad("'response' must be an object")
    else {
      fam <- r$family
      if (!.is_scalar_string(fam) || !fam %in% names(.family_params))
        bad("'response.family' must be one of ", paste(names(.family_params), collapse = ", "),
            if (.is_scalar_string(fam)) paste0(", not '", fam, "'") else "")
      if (!.is_scalar_string(r$name) || !nzchar(r$name))
        bad("'response.name' must be a non-empty string")
      if (.is_scalar_string(fam) && fam %in% names(.family_params)) {
        needed <- .family_params[[fam]]
        allowed <- c("family", "name", "round", needed)
        for (k in setdiff(names(r), allowed))
          unknown("'response.", k, "' is not used by the ", fam,
                  " family; it would be silently ignored")
        for (k in needed) if (is.null(r[[k]]))
          bad("'response.", k, "' is required for the ", fam, " family")
        if (!is.null(r$sigma) && (!.is_scalar_number(r$sigma) || r$sigma <= 0))
          bad("'response.sigma' must be greater than 0")
        if (!is.null(r$beta) && (!.is_scalar_number(r$beta) || r$beta <= 0))
          bad("'response.beta' must be greater than 0")
        if (!is.null(r$phi) && (!.is_scalar_number(r$phi) || r$phi <= 0))
          bad("'response.phi' must be greater than 0")
        if (!is.null(r$shift) && !.is_scalar_number(r$shift))
          bad("'response.shift' must be a number")
        if (!is.null(r$thresholds)) {
          th <- r$thresholds
          if (!is.numeric(th) || !length(th) || anyNA(th))
            bad("'response.thresholds' must be a non-empty numeric array")
          else if (length(th) > 1 && any(diff(th) <= 0))
            bad("'response.thresholds' must be strictly increasing")
        }
        if (!is.null(r$round)) {
          if (!.is_whole(r$round) || r$round < 0)
            bad("'response.round' must be a whole number of at least 0")
          else if (!fam %in% .rounding_families)
            unknown("'response.round' has no effect for the ", fam,
                    " family, whose outcome is already an integer")
        }
      }
    }
  }

  if (length(soft)) warning(paste0("in this design specification:\n  - ",
                                   paste(soft, collapse = "\n  - ")), call. = FALSE)
  if (length(problems))
    stop(paste0("invalid design specification:\n  - ", paste(problems, collapse = "\n  - ")),
         call. = FALSE)
  invisible(spec)
}
