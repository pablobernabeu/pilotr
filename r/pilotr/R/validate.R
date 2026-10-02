# Specification validation and version negotiation.
#
# load_spec() was a bare jsonlite::fromJSON call: no validation and no version check. A strict
# draft-07 schema has shipped at spec/design.schema.json since 0.1, and no code path consulted
# it, so a specification with a misspelled field loaded silently and kept the misspelling.
#
# That matters more than a typo usually would, because several of the ways a specification can
# be wrong produce plausible data and no error at all:
#
#   * A mistyped coefficient key resolves to no column, so that effect is silently set to zero.
#     A spec whose focal effect is spelled "cnod" for "cond" generates exactly the data of
#     a null design, and reports success.
#   * varies_by took anything other than "subject" as item-level, so a per-trial predictor was
#     silently given one value per item.
#   * A response parameter belonging to another family is silently ignored.
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
#
# Only well-formed parts are inspected. This scan runs before the per-field checks, so a
# malformed `response`, `predictors` or `random` value reached it first and stopped with a base
# error such as "$ operator is invalid for atomic vectors"; skipping it here leaves the field's
# own check to report it in the validator's words.
.spec_0_3_features <- function(spec) {
  found <- character(0)
  preds <- spec[["predictors"]]
  if (is.list(preds) && is.null(names(preds))) for (p in preds) {
    if (!is.list(p)) next
    if (identical(p[["varies_by"]], "observation"))
      found <- c(found, 'predictors varies_by "observation"')
    if (!is.null(p[["dist"]])) found <- c(found, "predictors dist")
    rel <- p[["reliability"]]
    if (!is.null(rel) && !(is.numeric(rel) && identical(as.numeric(rel), 1)))
      found <- c(found, "predictors reliability")
  }
  resp <- spec[["response"]]
  if (is.list(resp) && identical(resp[["family"]], "exgaussian"))
    found <- c(found, 'the "exgaussian" family')
  rs <- spec[["random"]]
  if (is.list(rs)) for (g in names(rs)) {
    re <- rs[[g]]
    if (!is.list(re)) next
    if (!is.null(re[["correlated"]])) found <- c(found, "random correlated")
    if (any(grepl(":", names(re[["slopes"]]), fixed = TRUE)))
      found <- c(found, "interaction random slopes")
  }
  unique(found)
}

# The largest seed both twins read exactly. jsonlite reads a JSON integer as a double, which holds
# every integer up to 2^53 - 1 and no longer every one beyond it, while Python reads an exact
# integer of any size. A larger seed was therefore two different numbers in the two readers, which
# then seeded different streams.
.MAX_EXACT_SEED <- 9007199254740991

# The first key that appears twice within one object of a parsed specification, or NULL when every
# key is unique. JSON leaves the meaning of a repeated key undefined (RFC 8259, section 4), and the
# two readers settle it differently: jsonlite keeps both entries, and a lookup by name finds the
# first of them however often it is made, while Python's json keeps only the last. A coefficient
# given twice was thus applied twice at its first value in R and once at its second in Python. The
# members of an object are searched before the object itself, which is the order in which Python's
# parser completes them, so that both twins name the same key.
.repeated_key <- function(x) {
  if (!is.list(x)) return(NULL)
  for (el in x) {
    k <- .repeated_key(el)
    if (!is.null(k)) return(k)
  }
  nms <- names(x)
  if (is.null(nms) || !anyDuplicated(nms)) return(NULL)
  nms[duplicated(nms)][[1L]]
}

.repeated_key_message <- function(key) {
  sprintf(paste0("the specification repeats the key '%s' within one object; JSON leaves a ",
                 "repeated key undefined, and R and Python read it differently"), key)
}

# Not `%||%`: base R gained that operator in 4.4.0, above the 4.0.0 pilotr declares, so defining
# it here would shadow the base version on new R and be the only definition on old R.
.orelse <- function(a, b) if (is.null(a)) b else a

# Resolve a specification argument to a validated list. Every public entry point calls this
# exactly once, so that a path is read once and validation runs once for a whole replicate loop.
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
# is refused, because a leftover parameter from another family is usually a half-finished edit
# and silently dropping it would hide the mistake.
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
#' data and no error at all. A mistyped coefficient key resolves to no column and so silently
#' sets that effect to zero, which generates exactly the data of a null design and reports
#' success. A response parameter left over from another family is ignored. Neither is
#' detectable in the output, which is why they are refused here.
#'
#' Version negotiation covers the other direction. A specification that uses a feature
#' introduced in 0.3 is read differently by a 0.2 implementation, so it must declare 0.3 or
#' later. A specification declaring a version newer than this implementation is refused
#' outright. A specification with no `spec_version` is treated as 0.2, which is
#' what every specification written before the field existed is.
#'
#' Every field is read by its exact name, as the 'Python' twin reads it, and a list that repeats
#' a name within one object is refused, as [load_spec()] refuses a file that does. A seed must
#' be a whole number within plus or minus 2^53 - 1, the range in which both twins read a JSON
#' integer exactly.
#'
#' @param spec A design specification (path or list).
#' @param strict Whether an unrecognised field is an error (the default) or a warning. Set
#'   `FALSE` to load a specification carrying private annotations, accepting that a misspelled
#'   field will then be ignored in silence.
#' @return The specification, invisibly, so that the call can be chained.
#' @examples
#' spec <- build_spec(list(name = "demo", seed = 1, design_kind = "between",
#'   factor_name = "group", lev1 = "a", lev2 = "b", n_subject = 20,
#'   intercept = 0, effect = 0.5, family = "gaussian", resp_name = "", sigma = 1))
#' validate_spec(spec)
#'
#' # A mistyped coefficient key is refused, where it used to pass as a zero effect.
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

  # Every field below is read with `[[ ]]`, never `$`. `$` matches a name partially, so it read
  # `random$subject` off an entry named `subject_site` and `units$item` off a misspelt `items`,
  # and the two twins then disagreed about what the specification said.

  # ---- repeated keys ----
  # load_spec() refuses these as it parses; this catches a hand-built list.
  key <- .repeated_key(spec)
  if (!is.null(key)) bad(.repeated_key_message(key))

  # ---- version ----
  declared <- .orelse(spec[["spec_version"]], "0.2")
  sv <- .parse_version(.SPEC_VERSION)
  # A version is one string or one number. Anything else used to reach .parse_version(), where an
  # empty array or object stopped with "subscript out of bounds" and a two-element array was read
  # by its first element, while the Python twin called both malformed.
  if (!.is_scalar_string(declared) && !.is_scalar_number(declared)) {
    bad("'spec_version' must be a single string of the form 'major.minor'")
  } else if (is.null(dv <- .parse_version(declared))) {
    bad("spec_version '", declared, "' is not of the form 'major.minor'")
  } else {
    # The version as pilotr read it, so that the two engines report the same thing about a JSON
    # number they render differently.
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

  if (!is.null(spec[["name"]]) && !.is_scalar_string(spec[["name"]]))
    bad("'name' must be a single string")
  seed <- spec[["seed"]]
  if (!is.null(seed) && !(.is_whole(seed) && abs(seed) <= .MAX_EXACT_SEED))
    bad("'seed' must be a whole number between -9007199254740991 and 9007199254740991 ",
        "(2^53 - 1), the range in which R and Python read a JSON integer exactly")

  # ---- units ----
  u <- spec[["units"]]
  has_item <- FALSE
  if (!is.null(u)) {
    if (!is.list(u) || is.null(names(u))) bad("'units' must be an object")
    else {
      for (k in setdiff(names(u), c("subject", "item")))
        unknown("unknown unit '", k, "'; only 'subject' and 'item' exist")
      if (is.null(u[["subject"]])) bad("'units.subject' is required")
      has_item <- !is.null(u[["item"]])
      for (nm in intersect(names(u), c("subject", "item"))) {
        un <- u[[nm]]
        # The list test here, where an object is elsewhere tested by names(): an empty JSON
        # object arrives as an unnamed empty list, and names() would call it the wrong shape
        # where the twin reports the missing n. Without this guard the read of `n` below threw a
        # base error such as "$ operator is invalid for atomic vectors" where Python said what
        # was wrong.
        if (!is.list(un)) { bad("'units.", nm, "' must be an object"); next }
        for (k in setdiff(names(un), c("n", "per_subject")))
          unknown("unknown field 'units.", nm, ".", k, "'")
        n <- un[["n"]]; per <- un[["per_subject"]]
        if (!.is_whole(n) || n < 1) bad("'units.", nm, ".n' must be a whole number of at least 1")
        if (!is.null(per)) {
          if (!identical(nm, "item"))
            bad("'per_subject' belongs to 'units.item', not 'units.", nm, "'")
          else if (!.is_whole(per) || per < 1)
            bad("'units.item.per_subject' must be a whole number of at least 1")
          else if (.is_whole(n) && per > n)
            bad("'units.item.per_subject' (", per, ") cannot exceed the number of items (", n, ")")
        }
      }
    }
  }

  # ---- factors ----
  contrast_cols <- character(0)
  factors <- spec[["factors"]]
  if (!is.null(factors)) {
    if (!is.list(factors) || !is.null(names(factors)))
      bad("'factors' must be an array of factor objects")
    else for (i in seq_along(factors)) {
      f <- factors[[i]]; where <- paste0("factors[", i, "]")
      if (!is.list(f) || is.null(names(f))) { bad(where, " must be an object"); next }
      for (k in setdiff(names(f), c("name", "levels", "contrasts", "vary_within", "between")))
        unknown("unknown field '", where, ".", k, "'")
      if (!.is_scalar_string(f[["name"]])) bad(where, ".name must be a single string")
      levels <- f[["levels"]]
      nlev <- length(levels)
      # A JSON null among the levels arrives as NA, which R would have used as a level label.
      if (!is.character(levels) || nlev < 2 || anyNA(levels))
        bad(where, ".levels must be an array of at least two strings")
      contrasts <- f[["contrasts"]]
      if (!is.list(contrasts) || is.null(names(contrasts)) || !length(contrasts))
        bad(where, ".contrasts must be a non-empty object mapping contrast columns to one value per level")
      else for (cn in names(contrasts)) {
        contrast_cols <- c(contrast_cols, cn)
        v <- contrasts[[cn]]
        if (!is.numeric(v) || anyNA(v)) bad(where, ".contrasts.", cn, " must be numeric")
        else if (nlev >= 2 && length(v) != nlev)
          bad(where, ".contrasts.", cn, " has ", length(v), " value(s) but the factor has ",
              nlev, " level(s)")
      }
      # A single string is accepted where an array belongs. pilotr's own spec_json() emitted that
      # form before 0.3, because a blanket auto_unbox collapsed every one-element array, so
      # refusing it would mean refusing files pilotr itself wrote. The reading is unambiguous and
      # both engines already treat the two alike.
      vw <- f[["vary_within"]]
      if (!is.null(vw)) {
        if (!is.character(vw) || !length(vw) || anyNA(vw))
          bad(where, ".vary_within must be a unit name or an array of unit names")
        else for (w in vw) {
          if (!w %in% c("subject", "item")) bad(where, ".vary_within contains '", w,
                                               "'; only 'subject' and 'item' are allowed")
          else if (identical(w, "item") && !has_item)
            bad(where, ".vary_within names 'item' but the design has no item unit")
        }
      }
      bt <- f[["between"]]
      if (!is.null(bt)) {
        if (!.is_scalar_string(bt) || !bt %in% c("subject", "item"))
          bad(where, ".between must be 'subject' or 'item'")
        else if (identical(bt, "item") && !has_item)
          bad(where, ".between is 'item' but the design has no item unit")
      }
      if (is.null(vw) && is.null(bt))
        bad(where, " must set either 'vary_within' or 'between'")
    }
  }

  # ---- predictors ----
  pred_names <- character(0)
  predictors <- spec[["predictors"]]
  if (!is.null(predictors)) {
    if (!is.list(predictors) || !is.null(names(predictors)))
      bad("'predictors' must be an array of predictor objects")
    else for (i in seq_along(predictors)) {
      p <- predictors[[i]]; where <- paste0("predictors[", i, "]")
      if (!is.list(p) || is.null(names(p))) { bad(where, " must be an object"); next }
      for (k in setdiff(names(p), c("name", "varies_by", "mean", "sd", "dist",
                                    "min", "max", "reliability")))
        unknown("unknown field '", where, ".", k, "'")
      if (!.is_scalar_string(p[["name"]])) bad(where, ".name must be a single string")
      else pred_names <- c(pred_names, p[["name"]])
      vb <- p[["varies_by"]]
      # Only a string is quoted back: pasting in an array split the one problem into several,
      # and the two twins render other values differently.
      if (!.is_scalar_string(vb) || !vb %in% c("subject", "item", "observation"))
        bad(where, ".varies_by must be 'subject', 'item' or 'observation'",
            if (.is_scalar_string(vb)) paste0(", not '", vb, "'") else "")
      else if (identical(vb, "item") && !has_item)
        bad(where, ".varies_by is 'item' but the design has no item unit")
      dist <- .orelse(p[["dist"]], "normal")
      if (!.is_scalar_string(dist) || !dist %in% c("normal", "uniform"))
        bad(where, ".dist must be 'normal' or 'uniform'")
      else if (identical(dist, "uniform")) {
        if (!.is_scalar_number(p[["min"]]) || !.is_scalar_number(p[["max"]]))
          bad(where, " uses dist 'uniform' and so needs numeric 'min' and 'max'")
        else if (p[["min"]] >= p[["max"]]) bad(where, ".min must be less than ", where, ".max")
        for (k in intersect(names(p), c("mean", "sd")))
          unknown(where, ".", k, " is ignored when dist is 'uniform'")
      } else {
        for (k in intersect(names(p), c("min", "max")))
          unknown(where, ".", k, " is ignored when dist is 'normal'")
        if (!is.null(p[["mean"]]) && !.is_scalar_number(p[["mean"]]))
          bad(where, ".mean must be a number")
        if (!is.null(p[["sd"]]) && (!.is_scalar_number(p[["sd"]]) || p[["sd"]] < 0))
          bad(where, ".sd must be a number of at least 0")
      }
      rel <- p[["reliability"]]
      if (!is.null(rel) && (!.is_scalar_number(rel) || rel <= 0 || rel > 1))
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
  fx <- spec[["fixed"]]
  if (!is.null(fx)) {
    if (!is.list(fx) || is.null(names(fx))) bad("'fixed' must be an object")
    else {
      for (k in setdiff(names(fx), c("intercept", "coefficients")))
        unknown("unknown field 'fixed.", k, "'")
      if (!.is_scalar_number(fx[["intercept"]])) bad("'fixed.intercept' must be a single number")
      coeffs <- fx[["coefficients"]]
      if (is.null(coeffs)) bad("'fixed.coefficients' is required (use {} for none)")
      else if (!is.list(coeffs)) bad("'fixed.coefficients' must be an object")
      else for (k in names(coeffs)) {
        if (!.is_scalar_number(coeffs[[k]]))
          bad("'fixed.coefficients.", k, "' must be a single number")
        check_key(k, "fixed.coefficients")
      }
    }
  }

  # ---- random ----
  rs <- spec[["random"]]
  if (!is.null(rs) && length(rs)) {
    if (!is.list(rs) || is.null(names(rs)))
      bad("'random' must be an object keyed by grouping factor")
    else for (g in names(rs)) {
      re <- rs[[g]]; where <- paste0("random.", g)
      if (!is.list(re) || is.null(names(re))) { bad(where, " must be an object"); next }
      for (k in setdiff(names(re), c("intercept_sd", "slopes", "correlations", "correlated",
                                     "over", "n")))
        unknown("unknown field '", where, ".", k, "'")
      if (!.is_scalar_number(re[["intercept_sd"]]) || re[["intercept_sd"]] < 0)
        bad(where, ".intercept_sd is required and must be at least 0")
      slopes <- re[["slopes"]]
      cols <- c("intercept", names(slopes))
      if (!is.null(slopes)) {
        if (!is.list(slopes)) bad(where, ".slopes must be an object")
        else for (k in names(slopes)) {
          if (!.is_scalar_number(slopes[[k]]) || slopes[[k]] < 0)
            bad(where, ".slopes.", k, " must be a number of at least 0")
          check_key(k, paste0(where, ".slopes"))
        }
      }
      cors <- re[["correlations"]]
      if (!is.null(cors)) {
        if (!is.list(cors)) bad(where, ".correlations must be an object")
        else for (k in names(cors)) {
          v <- cors[[k]]
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
      correlated <- re[["correlated"]]
      # `"correlated": [null]` arrives as NA, a logical that isTRUE() then read as FALSE.
      if (!is.null(correlated) &&
          !(is.logical(correlated) && length(correlated) == 1L && !is.na(correlated)))
        bad(where, ".correlated must be TRUE or FALSE")
      if (isTRUE(identical(correlated, FALSE)) && length(cors))
        bad(where, " sets correlated = false but also supplies correlations; one of the two has to go")
      if (g %in% c("subject", "item")) {
        for (k in intersect(names(re), c("over", "n")))
          bad(where, ".", k, " applies only to an extra grouping factor, not to '", g, "'")
      } else {
        over <- re[["over"]]
        if (!.is_scalar_string(over) || !over %in% c("subject", "item"))
          bad(where, ".over is required for an extra grouping factor and must be 'subject' or 'item'")
        else if (identical(over, "item") && !has_item)
          bad(where, ".over is 'item' but the design has no item unit")
        if (!.is_whole(re[["n"]]) || re[["n"]] < 1)
          bad(where, ".n is required for an extra grouping factor and must be a whole number of at least 1")
      }
    }
  }

  # ---- response ----
  r <- spec[["response"]]
  if (!is.null(r)) {
    if (!is.list(r) || is.null(names(r))) bad("'response' must be an object")
    else {
      fam <- r[["family"]]
      if (!.is_scalar_string(fam) || !fam %in% names(.family_params))
        bad("'response.family' must be one of ", paste(names(.family_params), collapse = ", "),
            if (.is_scalar_string(fam)) paste0(", not '", fam, "'") else "")
      if (!.is_scalar_string(r[["name"]]) || !nzchar(r[["name"]]))
        bad("'response.name' must be a non-empty string")
      if (.is_scalar_string(fam) && fam %in% names(.family_params)) {
        needed <- .family_params[[fam]]
        allowed <- c("family", "name", "round", needed)
        for (k in setdiff(names(r), allowed))
          unknown("'response.", k, "' is not used by the ", fam,
                  " family; it would be silently ignored")
        for (k in needed) if (is.null(r[[k]]))
          bad("'response.", k, "' is required for the ", fam, " family")
        for (k in c("sigma", "beta", "phi"))
          if (!is.null(r[[k]]) && (!.is_scalar_number(r[[k]]) || r[[k]] <= 0))
            bad("'response.", k, "' must be greater than 0")
        if (!is.null(r[["shift"]]) && !.is_scalar_number(r[["shift"]]))
          bad("'response.shift' must be a number")
        th <- r[["thresholds"]]
        if (!is.null(th)) {
          if (!is.numeric(th) || !length(th) || anyNA(th))
            bad("'response.thresholds' must be a non-empty numeric array")
          else if (length(th) > 1 && any(diff(th) <= 0))
            bad("'response.thresholds' must be strictly increasing")
        }
        if (!is.null(r[["round"]])) {
          if (!.is_whole(r[["round"]]) || r[["round"]] < 0)
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
