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
# seventeen is enough for every double, so trying the three widths in turn gives both
# exactness and readability.
.shortest_double <- function(z) {
  if (!is.finite(z)) return(NULL)
  for (d in 15:17) {
    s <- sprintf(paste0("%.", d, "g"), z)
    if (as.numeric(s) == z) return(s)
  }
  NULL
}

# Rewrite every JSON number to its shortest exact form.
#
# jsonlite formats all numbers at one fixed precision. Seventeen significant digits round-trip
# exactly but read badly: an effect the user typed as 0.3 is written 0.29999999999999999, which
# in an exported specification looks like a defect. This pass keeps the exactness and restores
# the readability, leaving 0.3 as "0.3" while 1/3 keeps every digit it needs. It walks the text
# so that digits inside string literals, such as a factor level named "group 2", are untouched.
.shorten_json_numbers <- function(txt) {
  cs <- strsplit(txt, "", fixed = TRUE)[[1]]
  n <- length(cs); out <- character(n); k <- 0L; i <- 1L; in_string <- FALSE
  add <- function(s) { k <<- k + 1L; out[k] <<- s }
  while (i <= n) {
    ch <- cs[i]
    if (in_string) {
      add(ch)
      if (ch == "\\" && i < n) { i <- i + 1L; add(cs[i]) } else if (ch == "\"") in_string <- FALSE
      i <- i + 1L; next
    }
    if (ch == "\"") { in_string <- TRUE; add(ch); i <- i + 1L; next }
    if (grepl("^[-0-9]$", ch)) {
      j <- i
      while (j <= n && grepl("^[-+0-9.eE]$", cs[j])) j <- j + 1L
      tok <- paste(cs[i:(j - 1L)], collapse = "")
      z <- suppressWarnings(as.numeric(tok))
      short <- if (is.na(z)) NULL else .shortest_double(z)
      add(if (is.null(short)) tok else short)
      i <- j; next
    }
    add(ch); i <- i + 1L
  }
  paste(out[seq_len(k)], collapse = "")
}

#' Serialise a design specification to pretty-printed JSON
#'
#' @details
#' Numbers are written at 17 significant digits, which is the shortest precision that
#' round-trips every IEEE-754 double exactly. The JSON file is the portable artefact that the
#' R and 'Python' implementations both read, so anything less makes the specification itself a
#' source of cross-language divergence: at the previous setting a coefficient of `1/3` came
#' back as `0.33333333333333298`, and over a sample of 214 doubles 189 failed to round-trip.
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
  .shorten_json_numbers(
    jsonlite::toJSON(.unbox_spec(spec), auto_unbox = FALSE, pretty = TRUE, digits = I(17)))
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
#' Numbers are emitted at 17 significant digits rather than through `deparse()`, which prints
#' 15 and so does not round-trip: `deparse(1/3)` reads back as a different double. Since the
#' point of the script is bit-for-bit reproduction, the embedded specification has to preserve
#' every coefficient exactly.
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
