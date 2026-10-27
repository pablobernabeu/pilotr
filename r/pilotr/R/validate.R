# Specification validation and version negotiation.
#
# By default, load_spec() and simulate_design() pass a specification through validate_spec(),
# which enforces the rules spec/design.schema.json states and the rules that relate one value to
# another. Without it, a specification with a misspelled field would load silently and keep the
# misspelling.
#
# That matters more than a typo usually would, because several of the ways a specification can
# be wrong produce plausible data and no error at all:
#
#   * A mistyped coefficient key resolves to no column, so that effect is silently set to zero.
#     A spec whose focal effect is spelled "cnod" for "cond" generates exactly the data of
#     a null design, and reports success.
#   * A response parameter belonging to another family is silently ignored.
#
# It also matters for version negotiation. A 0.3 specification opened by a 0.2 implementation,
# an un-upgraded Python twin, or a cached browser build produces different and wrong data while
# reporting success. A 0.2 implementation has no version check and cannot be fixed
# retrospectively, but a specification that uses a 0.3 feature can be made to say so, and every
# implementation from 0.3 onwards refuses what it does not understand.

# The newest specification version this implementation understands. What it writes declares the
# lowest version that reads it (.required_spec_version()), so that a specification needing nothing
# newer stays readable by the releases that know only 0.3.
.SPEC_VERSION <- "0.4"

# Features introduced in 0.3. A specification using any of them is read differently by a 0.2
# implementation, so it has to declare 0.3 or later.
#
# Only well-formed parts are inspected. This scan runs before the per-field checks, so a
# malformed `response`, `predictors` or `random` value reaches it first, and inspecting it would
# stop with a base error such as "subscript out of bounds". Skipping it here leaves the field's
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

# The names of the well-formed factors between each of `units`, kept for a unit only when it has
# two or more. Up to version 0.3, each between factor was assigned to blocks of its unit on its
# own, so two over one unit fell into the same or overlapping blocks. Version 0.4 allocates them
# jointly (simulate_design()). .spec_0_4_features() and .check_names() both read this, so the
# feature that needs 0.4 and the refusal below 0.4 cannot drift apart.
.shared_between <- function(factors, units) {
  out <- list()
  if (!is.list(factors) || !is.null(names(factors))) return(out)
  for (unit in units) {
    on_unit <- character(0)
    for (f in factors)
      if (is.list(f) && !is.null(names(f)) && is.null(f[["vary_within"]]) &&
          .is_name(f[["name"]]) && identical(f[["between"]], unit))
        on_unit <- c(on_unit, f[["name"]])
    if (length(on_unit) > 1L) out[[unit]] <- on_unit
  }
  out
}

# Features introduced in 0.4. A specification using any of them is read differently by an earlier
# implementation, so it has to declare 0.4 or later. One with several between factors on one unit
# that declares less is refused by .check_names(), in words that name both remedies, and one with
# a field new in 0.4 by validate_spec(). Only well-formed parts are inspected, as for the 0.3
# features.
.spec_0_4_features <- function(spec) {
  c(if (length(.shared_between(spec[["factors"]], c("subject", "item"))))
      "several between factors on one unit",
    .spec_0_4_fields(spec))
}

# The fields introduced in 0.4. A release that reads no version above 0.3 refuses such a field as
# unknown, or ignores it when validation is lenient, as releases before 0.3 always do. Ignored,
# `counterbalance` leaves every subject-item pair crossed with every level. It counts whatever its
# value, since `false` is as unknown to those releases as `true`.
.spec_0_4_fields <- function(spec) {
  factors <- spec[["factors"]]
  if (is.list(factors) && is.null(names(factors)))
    for (f in factors)
      if (is.list(f) && !is.null(names(f)) && !is.null(f[["counterbalance"]]))
        return("factors counterbalance")
  character(0)
}

# The lowest specification version that reads `spec` as written, never below 0.3. build_spec() and
# spec_from_model() declare it, so that what they write stays readable by the releases that read
# 0.3 (R 0.3.0 and 0.3.1, Python 0.3.0) unless it uses a 0.4 feature. A lower version would gain
# no reader, since releases before 0.3 ignore the field.
.required_spec_version <- function(spec) {
  if (length(.spec_0_4_features(spec))) "0.4" else "0.3"
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

# Resolve a specification argument to a validated list. The analysis functions, power_design(),
# power_mixed(), power_curve_mixed(), precision_design(), precision_curve(), sweep_spec(),
# response_variance(), calibrate_response(), brms_bridge() and generate_design_analysis(), call
# this once at their start, so that a path is read once and validation runs once for a whole
# replicate loop. The functions they hand the validated list to, simulate_design() (called with
# validate = FALSE), model_data() and model_formula(), read a path with load_spec() themselves
# and do not call this, so that the list is not validated again.
.as_spec <- function(spec, strict = TRUE) {
  if (is.character(spec)) return(load_spec(spec, validate = strict))
  validate_spec(spec, strict = strict)
  spec
}

.is_scalar_string <- function(x) is.character(x) && length(x) == 1L && !is.na(x)

# A name that reaches the data as a column name. A blank one passed .is_scalar_string(), and the
# spec then failed inside simulate_design() with base R's "replacement has length zero", which
# names neither the field nor the emptied control, where the twin wrote a column with no name.
.is_name <- function(x) .is_scalar_string(x) && nzchar(trimws(x))

# Joins names as "a", "a and b", "a, b and c".
.name_list <- function(x) {
  if (length(x) < 2L) return(paste0(x, collapse = ""))
  paste0(paste(x[-length(x)], collapse = ", "), " and ", x[length(x)])
}

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
#' Names are checked together as well as one by one. Two columns with one name, a contrast
#' column defined by two factors or named like another column, an interaction whose analysis
#' column already exists, a level listed twice, a correlation that pairs a term with itself or
#' gives one pair twice, a factor both between and within a unit, and a `random.item` entry in a
#' design without items are each refused, since each would otherwise validate and change the data
#' in silence. The rules are set out under Names in the specification. A within factor whose
#' `vary_within` leaves out a unit of the design draws a warning in either mode, since pilotr
#' crosses it with every unit anyway, and a specification declaring version 0.4 is refused for it.
#'
#' Several factors may be between one unit only in a specification declaring version 0.4 or
#' later, which allocates them jointly, one block of units per combination of their levels. Up to
#' version 0.3, pilotr assigned the levels of each between factor to blocks of its unit on its
#' own, so two between factors over one unit fell into the same or overlapping blocks. A 2 x 2
#' between-subjects design over 40 subjects gave cells of 20, 0, 0 and 20, and its effects could
#' not be estimated apart, so a specification declaring an earlier version is refused for them.
#' Every version reads the same design written as one between factor whose levels are the cells,
#' as the specification shows under Worked encodings, and [spec_from_model()] builds that factor
#' from a fitted pilot.
#'
#' From version 0.4, a within factor that sets `counterbalance` shows each subject each item at
#' one level, rotating the levels over subjects and items as a Latin square does. It needs an item
#' unit, and it is refused with `per_subject`, since every subject has to see every item. A
#' warning says when the subjects or the items are not a multiple of the number of lists, which is
#' the number of combinations of the counterbalanced factors' levels, so that the lists are
#' unequal. The simulated data name each subject's list in a column `list`, so no other column of
#' such a design may take that name.
#'
#' Version negotiation covers the other direction. A specification that uses a feature
#' introduced in 0.3 is read differently by a 0.2 implementation, so it must declare 0.3 or
#' later, and one that uses a feature introduced in 0.4 must declare 0.4. A specification
#' declaring a version newer than this implementation is refused outright. A specification with
#' no `spec_version` is treated as 0.2, which is what every specification written before the
#' field existed is.
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
#' # A mistyped coefficient key is refused rather than passed as a zero effect.
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
  # The declared version as a (major, minor) pair, kept for the rules that changed with 0.4. It
  # stays NULL when the version cannot be read, and those rules then apply as they did before 0.4.
  dv <- NULL
  # A version is one string or one number. Anything else is malformed, as the Python twin also
  # says, and is kept from .parse_version(), where an empty array or object would stop with
  # "subscript out of bounds" and a two-element array would be read by its first element.
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
    # Below 0.4 in the twin's tuple order, major version first.
    used <- .spec_0_4_fields(spec)
    if (length(used) && (dv[1] < 0L || (dv[1] == 0L && dv[2] < 4L)))
      bad("this specification uses ", paste(used, collapse = ", "),
          ", which requires spec_version \"0.4\", but declares ", shown, ". An implementation ",
          "that reads no version above 0.3 refuses such a field or ignores it, and ignoring it ",
          "generates different data.")
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
        # This tests for a list, where an object is elsewhere tested by names(). jsonlite reads an
        # empty JSON object as a named empty list, which the names() test also accepts, but an R
        # caller writes an empty object as list(), which has no names. The list test passes both
        # on to the report of the missing n, as the twin reports it for `{}`. It also passes an
        # empty JSON array, which jsonlite reads as list() and the twin as an empty object.
        # Without this guard, the read of `n` below threw a base error such as "$ operator is
        # invalid for atomic vectors" where Python said what was wrong.
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
  # Whether each subject sees a sample of the items, over which no factor can be counterbalanced.
  sampled <- has_item && is.list(u[["item"]]) && !is.null(u[["item"]][["per_subject"]])
  factors <- spec[["factors"]]
  if (!is.null(factors)) {
    if (!is.list(factors) || !is.null(names(factors)))
      bad("'factors' must be an array of factor objects")
    else for (i in seq_along(factors)) {
      f <- factors[[i]]; where <- paste0("factors[", i, "]")
      if (!is.list(f) || is.null(names(f))) { bad(where, " must be an object"); next }
      for (k in setdiff(names(f), c("name", "levels", "contrasts", "vary_within", "between",
                                    "counterbalance")))
        unknown("unknown field '", where, ".", k, "'")
      if (!.is_name(f[["name"]])) bad(where, ".name must be a non-empty string")
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
      # counterbalance (new in 0.4) shows each subject each item at one level, rotating the levels
      # over subjects and items, so it needs a within factor, an item unit and every item for
      # every subject. false means the same as leaving it out. `"counterbalance": [null]` arrives
      # as NA.
      cb <- f[["counterbalance"]]
      if (!is.null(cb)) {
        if (!is.logical(cb) || length(cb) != 1L || is.na(cb))
          bad(where, ".counterbalance must be true or false")
        else if (cb && is.null(vw))
          bad(where, ".counterbalance applies only to a within factor, one that sets vary_within")
        else if (cb && !has_item)
          bad(where, ".counterbalance rotates the levels over items, but the design has no item ",
              "unit")
        else if (cb && sampled)
          bad(where, ".counterbalance needs every subject to see every item, so it cannot be ",
              "combined with units.item.per_subject")
      }
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
      if (!.is_name(p[["name"]])) bad(where, ".name must be a non-empty string")
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
      if (!.is_name(g)) { bad("a 'random' grouping factor must have a non-empty name"); next }
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
        bad(where, ".correlated must be true or false")
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
      if (!.is_name(r[["name"]])) bad("'response.name' must be a non-empty string")
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
            bad("'response.thresholds' must be a number or a non-empty numeric array")
          else if (length(th) > 1 && any(diff(th) <= 0))
            bad("'response.thresholds' must be strictly increasing")
        }
        if (!is.null(r[["round"]])) {
          if (!.is_whole(r[["round"]]) || r[["round"]] < 0)
            bad("'response.round' must be a whole number of at least 0")
          else if (!fam %in% .rounding_families)
            unknown("'response.round' has no effect for the ", fam,
                    " family, whose outcome is already an integer")
          # Rounded to whole numbers, every proportion became 0 or 1, outside the family's
          # support. At one place or more, a proportion rounded onto 0 or 1 is moved inside the
          # interval (.beta_round() in core.R). At none, the interval holds no value to move to.
          else if (identical(fam, "beta") && r[["round"]] < 1)
            bad("'response.round' must be at least 1 for the beta family, because a proportion ",
                "rounded to whole numbers is 0 or 1, outside the family's support")
        }
      }
    }
  }

  # ---- names and structure ----
  found <- .check_names(spec, has_item, units_ok = is.list(u) && !is.null(names(u)), version = dv)
  problems <- c(problems, found$problems)
  soft <- c(soft, found$soft)

  if (length(soft)) warning(paste0("in this design specification:\n  - ",
                                   paste(soft, collapse = "\n  - ")), call. = FALSE)
  if (length(problems))
    stop(paste0("invalid design specification:\n  - ", paste(problems, collapse = "\n  - ")),
         call. = FALSE)
  invisible(spec)
}

# The rules that look at a specification's names together, run after the per-field checks and
# mirrored line for line by the twin's _check_names(). Each case they refuse would otherwise
# validate and then move an effect, rescale a variance or overwrite a column without a word. Where
# a column is overwritten, the twins also disagree, R replacing the earlier column and Python
# appending a second one under the same name, so one specification would give two tables.
# Only well-formed parts are inspected, since a malformed one is reported by its own check.
# `version` is the declared version as validate_spec() parsed it, or NULL when it could not be
# read. Returns the refusals and the warnings, in the order both twins report them.
.check_names <- function(spec, has_item, units_ok, version = NULL) {
  problems <- character(0); soft <- character(0)
  # Two rules changed with 0.4: several between factors on one unit became valid, and a
  # vary_within that leaves out a unit of the design became an error. The comparison is the
  # twin's tuple comparison, major first.
  from_0_4 <- !is.null(version) &&
    (version[1] > 0L || (version[1] == 0L && version[2] >= 4L))
  bad <- function(...) problems <<- c(problems, paste0(...))
  is_obj <- function(x) is.list(x) && !is.null(names(x))
  is_arr <- function(x) is.list(x) && is.null(names(x))
  factors <- if (is_arr(spec[["factors"]])) spec[["factors"]] else list()
  predictors <- if (is_arr(spec[["predictors"]])) spec[["predictors"]] else list()
  rs <- if (is_obj(spec[["random"]])) spec[["random"]] else list()
  groups <- Filter(.is_name, names(rs))
  r <- spec[["response"]]
  # Repeats in the order of their first appearance, as the twin lists them.
  repeats <- function(x) unique(x[x %in% x[duplicated(x)]])

  # A level listed twice put the declared effect into the data, since the simulator works by
  # position, while model_data(), which matches labels, gave every row the first level's value.
  for (i in seq_along(factors)) {
    f <- factors[[i]]
    if (!is_obj(f)) next
    lv <- f[["levels"]]
    if (is.character(lv) && !anyNA(lv))
      for (l in repeats(lv))
        bad("factors[", i, "].levels repeats '", l, "'; each level needs its own label")
    # With both fields each unit got one row per level, every one at the level that the between
    # assignment chose.
    if (!is.null(f[["vary_within"]]) && !is.null(f[["between"]]))
      bad("factors[", i, "] sets both 'vary_within' and 'between'; a factor has to set exactly ",
          "one of them")
  }

  # Up to version 0.3, the simulator assigned the levels of each between factor to blocks of its
  # unit on its own, so two such factors over one unit fell into the same or overlapping blocks.
  # A 2 x 2 between-subjects design over 40 subjects gave cells of 20, 0, 0 and 20, and neither
  # the second effect nor the interaction could then be estimated. Version 0.4 allocates such
  # factors jointly, so a specification has to declare it to have them. The other remedy, one
  # factor whose levels are the cells, keeps every cell in every version, and spec_from_model()
  # writes that encoding for a fitted pilot.
  shared <- if (from_0_4) list() else .shared_between(factors, c("subject", if (has_item) "item"))
  for (unit in names(shared)) {
    on_unit <- shared[[unit]]
    bad("the factors ", .name_list(sprintf("'%s'", on_unit)), " are ",
        if (length(on_unit) == 2L) "both" else "all", " between '", unit, "'. Below ",
        "specification version 0.4, pilotr assigns the levels of each between factor to blocks ",
        "of ", unit, "s on its own. The blocks of these factors then coincide or overlap, which ",
        "leaves some combinations of their levels without ", unit, "s and confounds their ",
        "effects. Declare spec_version \"0.4\", which allocates several between factors jointly. ",
        "Alternatively, encode the design as one between factor whose levels are the cells, give ",
        "it the contrast columns of these factors and key each interaction as 'a:b', as ",
        "spec_from_model() in R does for a fitted pilot. Earlier releases of pilotr read that ",
        "encoding too.")
  }

  # Every column of the simulated data, in the order the simulator writes them.
  cols <- character(0); roles <- character(0)
  claim <- function(nm, role) if (.is_name(nm)) {
    cols <<- c(cols, nm); roles <<- c(roles, role)
  }
  if (units_ok) {
    claim("subject", "the subject unit")
    if (has_item) claim("item", "the item unit")
    # From version 0.4, the data of a design with a counterbalanced factor name each subject's list.
    if (has_item && any(vapply(factors, function(f) is_obj(f) && isTRUE(f[["counterbalance"]]) &&
                                 !is.null(f[["vary_within"]]), logical(1))))
      claim("list", "the list column of a counterbalanced design")
  }
  for (g in setdiff(groups, c("subject", "item"))) claim(g, paste0("random.", g))
  for (i in seq_along(factors))
    if (is_obj(factors[[i]])) claim(factors[[i]][["name"]], paste0("factors[", i, "]"))
  for (i in seq_along(predictors))
    if (is_obj(predictors[[i]])) claim(predictors[[i]][["name"]], paste0("predictors[", i, "]"))
  if (is_obj(r)) claim(r[["name"]], "the response")
  for (nm in repeats(cols)) {
    who <- roles[cols == nm]
    # Two predictors alone are already reported as a duplicated predictor name.
    if (all(startsWith(who, "predictors["))) next
    bad("the name '", nm, "' is used for more than one column (", .name_list(who), "); every ",
        "unit, grouping factor, factor, predictor and the response needs a column of its own")
  }

  # Contrast columns are where the coefficients and slopes act. One defined by two factors
  # carried only the later factor's effect, and one named like a predictor took the
  # predictor's draws in place of the contrast. A contrast named after its own factor is
  # allowed: two-level specifications use it, and the factor's name is checked above.
  ccol <- character(0); cfac <- integer(0)
  for (i in seq_along(factors)) {
    cs <- if (is_obj(factors[[i]])) factors[[i]][["contrasts"]]
    if (is_obj(cs)) { ccol <- c(ccol, names(cs)); cfac <- c(cfac, rep(i, length(cs))) }
  }
  for (cn in unique(ccol)) {
    fs <- unique(cfac[ccol == cn])
    if (length(fs) > 1L)
      bad("the contrast column '", cn, "' is defined by more than one factor (",
          .name_list(sprintf("factors[%d]", fs)), "); each contrast column has to belong to a ",
          "single factor")
  }
  for (k in seq_along(ccol)) {
    cn <- ccol[k]; i <- cfac[k]
    if (identical(cn, factors[[i]][["name"]])) next
    who <- roles[cols == cn]
    if (length(who))
      bad("the contrast column '", cn, "' of factors[", i, "] is also the name of ",
          .name_list(who), "; a contrast column may share its own factor's name but no other ",
          "column's")
  }

  # model_data() writes an interaction key "a:b" to a product column "a_b", over any column
  # that already has that name.
  fx <- spec[["fixed"]]
  keys <- if (is_obj(fx) && is_obj(fx[["coefficients"]])) names(fx[["coefficients"]])
  for (g in groups) {
    re <- rs[[g]]
    if (is_obj(re) && is_obj(re[["slopes"]])) keys <- c(keys, names(re[["slopes"]]))
  }
  keys <- unique(keys[grepl(":", keys, fixed = TRUE)])
  for (j in seq_along(keys)) {
    col <- .us(keys[j])
    earlier <- keys[seq_len(j - 1L)]
    who <- c(roles[cols == col],
             sprintf("a contrast column of factors[%d]", unique(cfac[ccol == col])),
             sprintf("the analysis column of the interaction '%s'",
                     earlier[vapply(earlier, .us, character(1)) == col]))
    if (length(who))
      bad("the interaction '", keys[j], "' becomes the analysis column '", col,
          "', which is already the name of ", .name_list(who))
  }

  # A term paired with itself overwrote the unit diagonal of the correlation matrix, which
  # rescaled that term's variance, and a pair given twice kept whichever came later.
  for (g in groups) {
    re <- rs[[g]]
    cors <- if (is_obj(re)) re[["correlations"]]
    if (!is_obj(cors)) next
    seen <- character(0); seen_key <- character(0)
    for (k in names(cors)) {
      parts <- trimws(strsplit(gsub("~", ",", k), ",")[[1]])
      if (length(parts) != 2L) next
      if (identical(parts[1], parts[2])) {
        bad("random.", g, ".correlations key '", k, "' pairs '", parts[1], "' with itself; a ",
            "term's correlation with itself is always 1")
        next
      }
      # Both orders are looked up, which avoids sorting the pair, since R sorts by locale.
      ab <- paste(parts, collapse = "\r"); ba <- paste(rev(parts), collapse = "\r")
      hit <- match(c(ab, ba), seen)
      hit <- hit[!is.na(hit)]
      if (length(hit))
        bad("random.", g, ".correlations keys '", seen_key[hit[1]], "' and '", k, "' name the ",
            "same pair of terms; give each pair once")
      else { seen <- c(seen, ab); seen_key <- c(seen_key, k) }
    }
  }

  # The simulator drew no item effects without an item unit, while model_formula() still
  # emitted an item term for the analysis.
  if (units_ok && !has_item && !is.null(rs[["item"]]))
    bad("random.item describes an item unit the design does not have; add units.item or ",
        "remove random.item")

  # pilotr crosses a within factor with every unit of the design, whatever vary_within lists,
  # so an incomplete list has never changed the data. It draws a warning below version 0.4 and
  # is refused from 0.4.
  if (units_ok) {
    design_units <- c("subject", if (has_item) "item")
    for (i in seq_along(factors)) {
      f <- factors[[i]]
      if (!is_obj(f) || !is.null(f[["between"]])) next
      vw <- f[["vary_within"]]
      if (!is.character(vw) || !length(vw) || anyNA(vw) || !all(vw %in% design_units)) next
      miss <- setdiff(design_units, vw)
      if (!length(miss)) next
      text <- paste0(
        "factors[", i, "].vary_within lists ", .name_list(sprintf("'%s'", unique(vw))),
        " but not ", .name_list(sprintf("'%s'", miss)), ". pilotr crosses a within factor ",
        "with every unit of the design, so it varies within ",
        .name_list(sprintf("'%s'", miss)), " as well; list every unit, or make it between ",
        "'item' if items carry it.")
      if (from_0_4)
        bad(text, " This specification declares ", version[1], ".", version[2],
            ", which refuses it.")
      else soft <- c(soft, paste0(text, " From spec_version 0.4 this is an error."))
    }
  }

  # Counterbalanced factors rotate the C combinations of their levels over C lists, subject s
  # following list (s - 1) mod C. The lists balance only when the subjects and the items are both
  # multiples of C, and otherwise differ by one subject or one item. That draws a warning, not a
  # refusal, since the rotation simulates such a design as written, and a sweep over item counts
  # passes through it. Partial crossing is refused above, so it is not warned about here.
  un <- spec[["units"]]
  rot <- Filter(function(f)
    is_obj(f) && isTRUE(f[["counterbalance"]]) && !is.null(f[["vary_within"]]) &&
      is.null(f[["between"]]) && .is_name(f[["name"]]) && is.character(f[["levels"]]) &&
      length(f[["levels"]]) >= 2L && !anyNA(f[["levels"]]), factors)
  if (length(rot) && units_ok && has_item && is_obj(un[["subject"]]) && is_obj(un[["item"]]) &&
      is.null(un[["item"]][["per_subject"]])) {
    n_subject <- un[["subject"]][["n"]]
    n_item <- un[["item"]][["n"]]
    lists <- prod(vapply(rot, function(f) length(f[["levels"]]), integer(1)))
    if (.is_whole(n_subject) && n_subject >= 1 && .is_whole(n_item) && n_item >= 1) {
      reasons <- c(
        if (n_subject %% lists != 0) sprintf(paste0(
          " The %.0f subjects are not a multiple of %.0f, so the lists differ in size by one ",
          "subject, and an item is seen in some conditions by one subject more than in others."),
          n_subject, lists),
        if (n_item %% lists != 0) sprintf(paste0(
          " The %.0f items are not a multiple of %.0f, so each list shows some conditions on one ",
          "item more than others."), n_item, lists))
      names_rot <- .name_list(sprintf("'%s'", vapply(rot, function(f) f[["name"]], character(1))))
      lead <- if (length(rot) == 1L)
        sprintf("the counterbalanced factor %s rotates its %.0f levels over %.0f lists", names_rot,
                lists, lists)
      else sprintf(paste0("the counterbalanced factors %s rotate the %.0f combinations of their ",
                          "levels over %.0f lists"), names_rot, lists, lists)
      if (length(reasons))
        soft <- c(soft, paste0(lead, ", and the lists are unequal.", paste(reasons, collapse = "")))
    }
  }

  list(problems = problems, soft = soft)
}
