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

# A seed as a specification stores it. A whole seed inside R's integer range stays an integer, as
# it always has been, so existing JSON and generated scripts are unchanged. A larger whole seed
# stays a double, which as.integer() turned into NA with a coercion warning, although every seed
# up to 2^53 - 1 is valid and spec_json() writes it exactly. Anything else is returned untouched,
# for validate_spec() to refuse in its own words, where as.integer() truncated 1.5 to 1 in silence.
.as_seed <- function(x) {
  if (!.is_whole(x)) return(x)
  if (abs(x) < 2^31) as.integer(x) else as.numeric(x)
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
#' @details
#' A response of the Gaussian, lognormal or shifted lognormal family is rounded to four
#' decimals, and a `beta` response is left unrounded. Rounding moves every proportion near 0 or 1
#' onto the nearest value inside the interval at that many decimals. At four decimals, 9.8% of the
#' 40,000 rows of a design with `phi = 0.8` landed on 0.0001 or 0.9999, and rows at such a bound
#' dominate a Beta likelihood (see [brms_bridge()]). Unrounded, the data keep the distribution the
#' specification describes. They can then differ between R and 'Python' in their last bits where
#' the two do not share a maths library, as the cross-language vignette sets out.
#'
#' The specification declares the lowest `spec_version` that reads it, and nothing built here
#' needs a version above 0.3, so the releases of either package that read 0.3 read it too.
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
    name = p$name, seed = .as_seed(p$seed),
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
  # The lowest version that reads the specification, first, where spec_json() has always written
  # it. Stamping .SPEC_VERSION would make every specification the apps export unreadable by the
  # releases that read 0.3, although nothing built here needs 0.4.
  c(list(spec_version = .required_spec_version(spec)), spec)
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

# The shortest decimal string that `reader` reads back as exactly this double, or NULL. Fifteen
# significant digits is enough for most numbers a user types, and seventeen is enough for every
# double under a reader that rounds correctly, so trying the shorter widths first gives both
# exactness and readability. NULL means a value that is not finite, or, for reader = "r", one
# whose fifteen- and sixteen-digit forms some build of R might read as a different double.
#
# reader = "json" serves spec_json(). Its document is read by jsonlite in load_spec() and by
# Python's json, and both round correctly, so each candidate is checked with jsonlite itself.
# as.numeric() will not do as the check. On a build with long-double arithmetic, which includes
# R on x86 Windows and Linux, it rounds a decimal twice, first to long double and then to double.
# It reads "-8.46847779583186" as -0x1.0efdc5258p+3, for one, where jsonlite reads its
# neighbour, -0x1.0efdc52580001p+3. Checking with it put about one number in 5,000 to 7,000 a
# unit in the last place away from the value R held. Seventeen digits always round-trip under
# correct rounding, so that width needs no check.
#
# reader = "r" serves numbers written as R source, which R's parser reads with the routine
# as.numeric() uses. A candidate must also pass .r_reads_back(), and nothing at seventeen digits
# can, because those digits exceed 2^53.
.shortest_double <- function(z, reader = c("json", "r")) {
  reader <- match.arg(reader)
  if (!is.finite(z)) return(NULL)
  for (d in 15:16) {
    s <- sprintf(paste0("%.", d, "g"), z)
    if (.json_reads_back(s, z) && (reader == "json" || .r_reads_back(s, z))) return(s)
  }
  if (reader == "json") sprintf("%.17g", z) else NULL
}

# Whether jsonlite, the reader of load_spec(), reads the decimal `s` as exactly the double `z`.
# jsonlite parses a number with the C library's strtod(), which rounds correctly, as Python's
# float() does.
.json_reads_back <- function(s, z) isTRUE(jsonlite::fromJSON(s) == z)

# Whether R's parser reads the decimal `s` as exactly the double `z` on every build, given that a
# correctly rounded reader does. R reads a decimal with R_strtod5() (src/main/util.c), which
# accumulates the digits as an integer, divides or multiplies it once by a power of ten built by
# repeated squaring, and converts the result to double.
#
# Without long doubles, that is a single rounding, and so a correct one, provided the integer
# stays within 2^53 and the power of ten within 10^22, the largest that is exact. With long
# doubles, the result is first rounded to the long double's significand, 64 bits on x86 and 113
# where long double has quadruple precision, and then to the double's 53. The second rounding
# goes astray only when the first lands exactly halfway between two doubles. The exact result
# is therefore worked out here, and the decimal is refused when it lies within one unit of a
# 64-bit significand of the midpoint. That is twice as far as a 64-bit rounding can move it.
#
# A decimal with k digits after the point lies a multiple of 1/5^k of the spacing of doubles
# from the double below it. Since 5^k is odd, it lies at least 1/(2 * 5^k) of that spacing from
# a midpoint, which for k = 4 is 1.6 units of a 64-bit significand. A value typed with up to
# four decimal places, such as 0.05 or 1.01, is therefore never refused. Of the 15- and 16-digit
# decimals that the other conditions admit, the margin refuses about one in 1,400. The check
# against as.numeric() confirms the outcome on the build writing the script.
.r_reads_back <- function(s, z) {
  if (!isTRUE(as.numeric(s) == z)) return(FALSE)
  if (z == 0) return(TRUE)
  t <- sub("^-", "", s)
  expn <- 0L
  if (grepl("e", t, fixed = TRUE)) {
    expn <- as.integer(sub("^.*e", "", t))
    t <- sub("e.*$", "", t)
  }
  if (grepl(".", t, fixed = TRUE)) expn <- expn - nchar(sub("^.*[.]", "", t))
  digits <- sub("^0+", "", sub(".", "", t, fixed = TRUE))
  # At most 2^53 = 9007199254740992, compared in two halves that every build reads exactly.
  n <- nchar(digits)
  if (n > 16L) return(FALSE)
  if (n == 16L) {
    hi <- as.numeric(substr(digits, 1L, 8L))
    if (hi > 90071992 || (hi == 90071992 && as.numeric(substr(digits, 9L, 16L)) > 54740992))
      return(FALSE)
  }
  if (abs(expn) > 22L) return(FALSE)

  m <- as.numeric(digits)
  p10 <- 1
  for (i in seq_len(abs(expn))) p10 <- p10 * 10
  a <- abs(z)
  # r = q - a, where q is the exact quotient or product, to far more than 64 bits.
  r <- if (expn < 0L) {
    ap <- .two_prod(a, p10)
    ((m - ap[1L]) - ap[2L]) / p10
  } else {
    mp <- .two_prod(m, p10)
    (mp[1L] - a) + mp[2L]
  }
  # q lies in [2^e, 2^(e + 1)), where doubles are 2^(e - 52) apart and the units of a 64-bit
  # significand 2^(e - 63). That is a's range, except when q lies just below a power of two.
  e <- floor(log2(a))
  if (2^e > a) e <- e - 1
  if (2^(e + 1) <= a) e <- e + 1
  if (r < 0 && a == 2^e) e <- e - 1
  # The midpoint on q's side lies half a spacing from a. A 64-bit rounding reaches it from within
  # half a unit, 2^(e - 64). The margin is twice that, and the error in r is far smaller still.
  abs(abs(r) - 2^(e - 53)) > 2^(e - 63)
}

# The product a * b as an unevaluated sum of two doubles, exactly: Dekker's (1971) algorithm,
# as Ogita, Rump and Oishi (2005) write it. Veltkamp's split, through the factor 2^27 + 1, cuts
# each factor into two halves whose pairwise products are exact.
.two_prod <- function(a, b) {
  p <- a * b
  ca <- 134217729 * a
  a1 <- ca - (ca - a)
  a2 <- a - a1
  cb <- 134217729 * b
  b1 <- cb - (cb - b)
  b2 <- b - b1
  c(p, a2 * b2 - (((p - a1 * b1) - a2 * b1) - a1 * b2))
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
#' A width is accepted only when 'jsonlite' reads it back as exactly the double it was written
#' from. [load_spec()] reads the file with 'jsonlite', which parses numbers with the C library's
#' `strtod()`, and the 'Python' package reads it with `json`. Both round correctly, so the two
#' packages recover the same double. Version 0.3.1 checked each width with R's `as.numeric()`,
#' which on a build with long-double arithmetic, such as R on x86 Windows and Linux, rounds a
#' decimal twice, first to long double and then to double. It accepted some 15- and 16-digit
#' forms that the readers of the file put one unit in the last place away.
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

# A double as R source that evaluates to the same bit pattern on every build. deparse() prints
# 15 significant digits, so deparse(1/3) reads back as a different double. A plain decimal is
# written where R's parser reads one exactly on every build, as it does any value typed with up
# to four decimal places. Any other number, a fitted value needing 17 digits for instance, is
# written as a call that has jsonlite read its decimal, since jsonlite rounds correctly on every
# build.
.num_literal <- function(z) {
  if (is.na(z)) return("NA_real_")
  if (!is.finite(z)) return(if (z > 0) "Inf" else "-Inf")
  short <- .shortest_double(z, "r")
  if (!is.null(short)) return(short)
  sprintf("jsonlite::fromJSON(\"%s\")", .shortest_double(z))
}

# Text as an R raw string, which the parser takes character for character, quotation marks and
# backslashes included (R >= 4.0.0). Dashes are added to the delimiter until its closing form, a
# parenthesis, the dashes and a quotation mark, occurs nowhere in the text.
.r_raw_string <- function(txt) {
  dashes <- ""
  while (grepl(paste0(")", dashes, "\""), txt, fixed = TRUE)) dashes <- paste0(dashes, "-")
  paste0("r\"", dashes, "(\n", txt, "\n)", dashes, "\"")
}

# Backtick-quote a list name unless it is already a syntactic R name. Coefficient and
# correlation keys such as "cond:z_freq" and "intercept,cond" are not.
.r_name <- function(n) {
  if (grepl("^[A-Za-z.][A-Za-z0-9._]*$", n) && !grepl("^\\.[0-9]", n)) n else paste0("`", n, "`")
}

# An R value as an R expression, built term by term so that every number goes through
# .num_literal(). deparse() would be shorter but silently rounds the numbers.
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
#' Embed the specification as the JSON text that [spec_json()] writes, so that the returned
#' script reproduces the design without any external file. This turns a design built in the
#' no-code application into a reproducible script. The application's Verify button runs that
#' script in a clean R session and checks that its data are identical, bit for bit, to the data
#' the application's current specification produces.
#'
#' @details
#' The script ends with a commented suggestion of the analysis that fits the design. A two-group
#' Gaussian design with one row per subject is sent to [power_design()], whose t-test takes no
#' other. An `ordinal` or `beta` outcome is sent to [generate_design_analysis()], which writes a
#' Bayesian design analysis on the link scale to `design_analysis.R`, since [power_mixed()] fits
#' such an outcome only on the response scale. Any other design is sent to [power_mixed()], and
#' only then does the installation line name `lme4` and `lmerTest`, which it needs.
#'
#' Since the point of the script is bit-for-bit reproduction, the embedded specification has to
#' preserve every coefficient exactly. The JSON text sits in a raw string that [load_spec()]
#' reads with 'jsonlite', which parses numbers with the C library's `strtod()`. That rounds
#' correctly, so each number reaches the script as the double it was on every build of R, and
#' the numbers stay as readable as in the specification file. Written as R literals, they would
#' be left to R's own parser. On a build without long-double arithmetic, such as R on Apple
#' silicon, that parser can read a 17-digit decimal one unit in the last place away. `deparse()`
#' is no alternative either, since it prints 15 significant digits and `deparse(1/3)` reads back
#' as a different double.
#'
#' @param spec A design specification (list), as produced by [build_spec()], or the path of a
#'   specification file for [load_spec()] to read.
#' @return A length-one character string containing a runnable R script that loads `pilotr`,
#'   embeds the specification and simulates the data, with the installation line and the
#'   suggested analysis as comments.
#' @examples
#' spec <- build_spec(list(name = "demo", seed = 1, design_kind = "between",
#'   factor_name = "group", lev1 = "a", lev2 = "b", n_subject = 20,
#'   intercept = 0, effect = 0.5, family = "gaussian", resp_name = "", sigma = 1))
#' cat(generate_r_script(spec))
#' @export
generate_r_script <- function(spec) {
  # A path used to reach spec_json() as a plain string, so the script embedded the path itself
  # and ran only where that file was. The file is now read, and the closing suggestion is chosen
  # from the design it holds.
  if (is.character(spec) && length(spec) == 1L) spec <- load_spec(spec)
  if (!is.list(spec))
    stop("`spec` must be a specification list or the path of a specification file",
         call. = FALSE)
  hint <- .script_hint(spec)
  paste0(
    "# Reproducible simulation exported by pilotr.\n",
    hint[["install"]],
    "library(pilotr)\n\n",
    "# The specification as spec_json() writes it. load_spec() reads it with jsonlite, which\n",
    "# recovers every number exactly on every build of R.\n",
    "spec <- load_spec(", .r_raw_string(spec_json(spec)), ")\n\n",
    "data <- simulate_design(spec)              # analysis-ready data frame\n",
    "# write.csv(data, \"data.csv\", row.names = FALSE)\n",
    hint[["analysis"]]
  )
}

# The installation line and the closing suggestion of an exported script, both chosen by the
# design. The script used to suggest power_mixed() for every design, beside an installation line
# that named none of the packages power_mixed() needs. power_design() takes the two-group Gaussian
# design with one row per subject, by its own rule. power_mixed() fits an ordinal or Beta outcome
# only on the response scale, while the specification's coefficients are on the logit scale.
# Those outcomes therefore get the Bayesian design analysis, written to a file with the first
# coefficient as the focal effect, as both apps' power tabs suggest. Called without a file, it
# returns the analysis as a string, and a script run with Rscript or source() discards that.
# Every other design gets power_mixed(), and the installation line then names lme4 and lmerTest.
.script_hint <- function(spec) {
  cran <- "# install.packages(\"pilotr\")   # from CRAN\n"
  if (is.null(.two_group_refusal(spec)))
    return(list(
      install = cran,
      analysis = "# pow <- power_design(spec, n_sims = 1000)   # power with Type S and Type M\n"))
  if (isTRUE(spec[["response"]][["family"]] %in% c("ordinal", "beta"))) {
    focal <- names(spec[["fixed"]][["coefficients"]])
    focal <- if (length(focal)) encodeString(focal[1L], quote = "\"") else "..."
    return(list(install = cran, analysis = sprintf(paste0(
      "# generate_design_analysis(spec, focal = %s, file = \"design_analysis.R\")   ",
      "# Bayesian design analysis on the link scale\n"), focal)))
  }
  list(install = paste0("# install.packages(c(\"pilotr\", \"lme4\", \"lmerTest\"))   # from CRAN; ",
                        "lme4 and lmerTest only for power_mixed()\n"),
       analysis = "# pow <- power_mixed(spec, n_sims = 200)   # needs lme4 and lmerTest\n")
}
