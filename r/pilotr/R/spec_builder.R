# Logic that converts GUI inputs into a portable spec, provided as package functions so that
# the application remains a thin client and this logic can be unit-tested. The code is plain R
# and does not depend on Shiny.

#' Default response-column name for a family
#'
#' @param family A response-family name, one of `"gaussian"`, `"lognormal"`,
#'   `"shifted_lognormal"`, `"bernoulli"`, `"poisson"`, `"ordinal"`, or
#'   `"beta"`.
#' @return The conventional response-column name for that family (for example `"RT"` for
#'   `"lognormal"` and `"shifted_lognormal"`), or `"outcome"` for an unrecognised family.
#' @examples
#' default_response_name("bernoulli")
#' @export
default_response_name <- function(family) {
  switch(family,
         gaussian = "score", lognormal = "RT", shifted_lognormal = "RT",
         bernoulli = "accuracy", poisson = "count", ordinal = "rating",
         beta = "proportion", "outcome")
}

#' Build a design specification from a flat list of design inputs
#'
#' Assemble a portable design specification (a plain list, serialisable with
#' [spec_json()]) from the flat set of inputs collected by the no-code
#' application: sample sizes, the two-level factor and its levels, the fixed
#' intercept and effect, the random-effect standard deviations for
#' within-subject and crossed designs, and the response family with its
#' parameters.
#'
#' @param p A named list of design inputs. Common fields are `name`, `seed`,
#'   `n_subject`, `design_kind` (`"between"` or `"within"`), `include_items`,
#'   `n_item`, `factor_name`, `lev1`, `lev2`, `intercept`, `effect`,
#'   `family`, `resp_name`, and family parameters such as `sigma`, `shift`,
#'   `thresholds`, or `phi`; within-design random effects use `subj_int_sd`,
#'   `subj_slope_sd`, `subj_corr`, `item_int_sd`, `item_slope_sd`, and
#'   `item_corr`.
#' @return A design specification as a nested list, ready for
#'   [simulate_design()], [spec_json()], or the
#'   power and precision functions.
#' @examples
#' build_spec(list(name = "demo", seed = 1, design_kind = "between",
#'   factor_name = "group", lev1 = "control", lev2 = "treatment", n_subject = 40,
#'   intercept = 100, effect = 5, family = "gaussian", resp_name = "", sigma = 10))
#' @export
build_spec <- function(p) {
  resp_name <- if (is.null(p$resp_name) || !nzchar(p$resp_name)) default_response_name(p$family) else p$resp_name

  factor <- list(name = p$factor_name, levels = c(p$lev1, p$lev2),
                 contrasts = list(effect = c(-0.5, 0.5)))
  if (identical(p$design_kind, "within")) {
    factor$vary_within <- if (isTRUE(p$include_items)) c("subject", "item") else c("subject")
  } else {
    factor$between <- "subject"
  }

  units <- list(subject = list(n = as.integer(p$n_subject)))
  if (identical(p$design_kind, "within") && isTRUE(p$include_items)) {
    units$item <- list(n = as.integer(p$n_item))
  }

  spec <- list(
    spec_version = .SPEC_VERSION,
    name = p$name, seed = as.integer(p$seed),
    units = units,
    factors = list(factor),
    fixed = list(intercept = p$intercept, coefficients = list(effect = p$effect))
  )

  if (identical(p$design_kind, "within")) {
    subj <- list(intercept_sd = p$subj_int_sd)
    if (isTRUE(p$subj_slope_sd > 0)) {
      subj$slopes <- list(effect = p$subj_slope_sd)
      subj$correlations <- list(`intercept,effect` = p$subj_corr)
    }
    spec$random <- list(subject = subj)
    if (isTRUE(p$include_items)) {
      item <- list(intercept_sd = p$item_int_sd)
      if (isTRUE(p$item_slope_sd > 0)) {
        item$slopes <- list(effect = p$item_slope_sd)
        item$correlations <- list(`intercept,effect` = p$item_corr)
      }
      spec$random$item <- item
    }
  }

  resp <- list(family = p$family, name = resp_name)
  if (p$family %in% c("gaussian", "lognormal", "shifted_lognormal")) { resp$sigma <- p$sigma; resp$round <- 4L }
  if (p$family == "shifted_lognormal") resp$shift <- p$shift
  if (p$family == "ordinal") resp$thresholds <- as.numeric(strsplit(gsub("\\s", "", p$thresholds), ",")[[1]])
  if (p$family == "beta") resp$phi <- p$phi
  spec$response <- resp
  spec
}

# Fields that stay JSON arrays even when they hold a single value. Everything else of length
# one is a JSON scalar. A blanket `auto_unbox = TRUE` could not draw this distinction: it
# collapsed a one-element `vary_within`, a single-level `levels`, a one-cut-point `thresholds`
# or a single-level contrast into a bare value, which violates design.schema.json and does not
# round-trip through load_spec().
.spec_array_fields <- c("levels", "vary_within", "thresholds")

# Fields that are JSON objects, so that an empty one emits `{}` rather than `[]`. R cannot tell
# an empty named list from an empty unnamed one, and `"random": []` fails the schema.
.spec_object_fields <- c("units", "subject", "item", "contrasts", "fixed", "coefficients",
                         "random", "slopes", "correlations", "response")

# Recursive walk marking genuine scalars for jsonlite::unbox(). `parent` carries the enclosing
# key so that the members of `contrasts` (one numeric array per contrast column) stay arrays.
.unbox_spec <- function(x, key = NULL, parent = NULL) {
  if (is.list(x)) {
    if (length(x) == 0L)
      return(if (!is.null(key) && key %in% .spec_object_fields)
        structure(list(), names = character(0)) else x)
    nms <- names(x)
    if (is.null(nms)) return(lapply(x, .unbox_spec, key = key, parent = parent))
    out <- lapply(seq_along(x), function(i) .unbox_spec(x[[i]], key = nms[i], parent = key))
    names(out) <- nms
    return(out)
  }
  keep_array <- identical(parent, "contrasts") ||
    (!is.null(key) && key %in% .spec_array_fields)
  if (keep_array || length(x) != 1) return(x)
  jsonlite::unbox(x)
}

# The shortest decimal string that reads back as exactly this double, or NULL for a
# non-finite value. Fifteen significant digits is enough for most numbers a user types, and
# seventeen is enough for every double, so trying the shorter widths first and falling back to
# seventeen gives both exactness and readability.
#
# The seventeen-digit form is returned without being checked, because on a build of R without
# long-double arithmetic the check would reject the one form that is certainly right. There
# `as.numeric()` accumulates the mantissa in a double before applying the decimal exponent, and
# seventeen digits overflow the 53-bit mantissa on the way, so "0.33333333333333331" reads back
# one unit in the last place below 1/3. Fifteen and sixteen digits stay within the mantissa and
# are correctly rounded everywhere, which is also why they are tried first: the shorter form is
# both the readable one and the one that survives a reader built without long doubles.
.shortest_double <- function(z) {
  if (!is.finite(z)) return(NULL)
  for (d in 15:16) {
    s <- sprintf(paste0("%.", d, "g"), z)
    if (as.numeric(s) == z) return(s)
  }
  sprintf("%.17g", z)
}

# The wrapper that carries an already-formatted number through jsonlite as a string.
#
# It is lengthened until it appears nowhere in the specification's own text, because a factor
# level or a label written to imitate it would otherwise be unwrapped into a number. Only '@'
# is ever added, so the tag never acquires a regular-expression metacharacter and
# .untag_json_numbers() can go on matching it literally.
.free_number_tag <- function(spec) {
  flat <- unlist(spec)
  text <- c(names(flat), as.character(flat))
  tag <- "@pilotr-number@"
  while (any(grepl(tag, text, fixed = TRUE))) tag <- paste0("@", tag)
  tag
}

# Format every finite double at its shortest exact width and wrap it in the tag, so that the
# number reaches the document as pilotr wrote it.
#
# jsonlite formats all numbers at one fixed precision. Seventeen significant digits round-trip
# exactly but read badly: an effect the user typed as 0.3 would be written 0.29999999999999999,
# which in an exported specification looks like a defect. Formatting each double here, rather
# than shortening jsonlite's output afterwards, keeps the readability without ever reading a
# number back from the text it was just written to. That read-back was wrong on a build of R
# without long-double arithmetic, where `as.numeric("0.33333333333333331")` is one unit in the
# last place below 1/3: the shortened form then recorded that wrong value, and the
# specification no longer held the coefficient the user set.
#
# Integers are left to jsonlite, which writes them exactly. So is any vector holding NA, NaN or
# an infinity, which has no decimal form to shorten and which jsonlite already renders.
.tag_numbers <- function(x, tag) {
  if (is.list(x)) return(lapply(x, .tag_numbers, tag = tag))
  if (!is.double(x) || !all(is.finite(x))) return(x)
  paste0(tag, vapply(x, .shortest_double, character(1)), tag)
}

# Unwrap the tagged numbers, dropping the quotation marks jsonlite put round them.
#
# The tag is already known not to occur in the specification's own text. Requiring the content
# between the two tags to be made of the characters a decimal number is written with is the
# second guard, and it is what keeps a malformed document from being unwrapped into one that
# no longer parses.
.untag_json_numbers <- function(txt, tag) {
  gsub(paste0("\"", tag, "([-+0-9.eE]+)", tag, "\""), "\\1", txt)
}

#' Serialise a design specification to pretty-printed JSON
#'
#' @details
#' Each number is written at the shortest width that reads back as exactly the same
#' IEEE-754 double, which is 15 significant digits for most values a user types and at most
#' 17 for the rest. The JSON file is the portable artefact that the R and 'Python'
#' implementations both read, so a fixed lower precision makes the specification itself a
#' source of cross-language divergence: at 15 digits throughout, a coefficient of `1/3` came
#' back as `0.33333333333333298`, and over a sample of 214 doubles 189 failed to round-trip.
#' Writing every number at 17 instead would round-trip but read badly, turning an effect
#' typed as 0.3 into `0.29999999999999999`.
#'
#' Preferring the shorter width is not only a matter of appearance. A build of R without
#' long-double arithmetic reads 17 significant digits inexactly, because the digits overflow
#' the 53-bit mantissa before the decimal exponent is applied, so the few values that need
#' that width are the ones such a build cannot recover. Anything shorter it reads correctly.
#'
#' @param spec A design specification (list), as produced by [build_spec()].
#' @return A length-one character string containing the specification as pretty-printed JSON,
#'   the portable artefact that the R and 'Python' packages both consume.
#' @examples
#' spec <- build_spec(list(name = "demo", seed = 1, design_kind = "between",
#'   factor_name = "group", lev1 = "a", lev2 = "b", n_subject = 20,
#'   intercept = 0, effect = 0.5, family = "gaussian", resp_name = "", sigma = 1))
#' cat(spec_json(spec))
#' @export
spec_json <- function(spec) {
  tag <- .free_number_tag(spec)
  .untag_json_numbers(
    jsonlite::toJSON(.unbox_spec(.tag_numbers(spec, tag)), auto_unbox = FALSE, pretty = TRUE,
                     digits = I(17)),
    tag)
}

# A double as an R source literal that parses back to the same bit pattern. deparse() prints
# 15 significant digits, so deparse(1/3) reads back as a different double.
.num_literal <- function(z) {
  if (is.na(z)) return("NA_real_")
  short <- .shortest_double(z)
  if (is.null(short)) return(if (z > 0) "Inf" else "-Inf")
  short
}

# Backtick-quote a list name unless it is already a syntactic R name. Coefficient and
# correlation keys such as "cond:z_freq" and "intercept,cond" are not.
.r_name <- function(n) {
  if (grepl("^[A-Za-z.][A-Za-z0-9._]*$", n) && !grepl("^\\.[0-9]", n)) n else paste0("`", n, "`")
}

# The specification as an R expression, built term by term so that every number keeps full
# precision. deparse() on the whole list would be shorter but silently rounds the numbers.
.r_literal <- function(x) {
  if (is.null(x)) return("NULL")
  if (is.list(x)) {
    nms <- names(x)
    parts <- vapply(seq_along(x), function(i) {
      v <- .r_literal(x[[i]])
      if (!is.null(nms) && nzchar(nms[i])) paste0(.r_name(nms[i]), " = ", v) else v
    }, character(1))
    return(paste0("list(", paste(parts, collapse = ", "), ")"))
  }
  v <- if (is.character(x)) encodeString(x, quote = "\"")
  else if (is.logical(x)) ifelse(is.na(x), "NA", ifelse(x, "TRUE", "FALSE"))
  else if (is.integer(x)) ifelse(is.na(x), "NA_integer_", paste0(format(x), "L"))
  else vapply(x, .num_literal, character(1))
  if (length(v) == 1) v else paste0("c(", paste(v, collapse = ", "), ")")
}

#' Generate a self-contained, reproducible R script from a specification
#'
#' Embed the specification as an R list literal, so that the returned script reproduces the
#' design without any external file. This turns a design built in the no-code application into
#' a reproducible script; the application's Verify button runs that script in a clean R session
#' and confirms that it reproduces the data bit-for-bit.
#'
#' @details
#' Numbers are emitted at the shortest width that reads back as the same double, rather than
#' through `deparse()`, which prints 15 significant digits and so does not round-trip:
#' `deparse(1/3)` reads back as a different double. Since the point of the script is
#' bit-for-bit reproduction, the embedded specification has to preserve every coefficient
#' exactly.
#'
#' @param spec A design specification (list), as produced by [build_spec()].
#' @return A length-one character string containing a runnable R script that loads `pilotr`,
#'   embeds the specification, and simulates the data.
#' @examples
#' spec <- build_spec(list(name = "demo", seed = 1, design_kind = "between",
#'   factor_name = "group", lev1 = "a", lev2 = "b", n_subject = 20,
#'   intercept = 0, effect = 0.5, family = "gaussian", resp_name = "", sigma = 1))
#' cat(generate_r_script(spec))
#' @export
generate_r_script <- function(spec) {
  paste0(
    "# Reproducible simulation exported by pilotr.\n",
    "# install.packages(\"pilotr\")   # once available; then run this script as-is.\n",
    "library(pilotr)\n\n",
    "spec <- ", .r_literal(spec), "\n\n",
    "data <- simulate_design(spec)              # analysis-ready data frame\n",
    "# write.csv(data, \"data.csv\", row.names = FALSE)\n",
    "# pow  <- power_mixed(spec, n_sims = 200)   # simulation-based power + Type S/M\n"
  )
}
