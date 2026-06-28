# Logic that converts GUI inputs into a portable spec, provided as package functions so that
# the application remains a thin client and this logic can be unit-tested. The code is plain R
# and does not depend on Shiny.

#' Default response-column name for a family
#'
#' @param family A response-family name, one of `"gaussian"`, `"shifted_lognormal"`,
#'   `"bernoulli"`, `"poisson"`, `"ordinal"`, or `"beta"`.
#' @return The conventional response-column name for that family (for example `"RT"` for
#'   `"shifted_lognormal"`), or `"outcome"` for an unrecognised family.
#' @examples
#' default_response_name("bernoulli")
#' @export
default_response_name <- function(family) {
  switch(family,
         gaussian = "score", shifted_lognormal = "RT",
         bernoulli = "accuracy", poisson = "count", ordinal = "rating",
         beta = "proportion", "outcome")
}

#' Build a design specification from a flat list of design inputs
#'
#' Assemble a portable design specification (a plain list, serialisable with [spec_json()])
#' from the flat set of inputs collected by the no-code application: sample sizes, the
#' two-level factor and its levels, the fixed intercept and effect, the random-effect
#' standard deviations for within-subject and crossed designs, and the response family with
#' its parameters.
#'
#' @param p A named list of design inputs. Common fields are `name`, `seed`, `n_subject`,
#'   `design_kind` (`"between"` or `"within"`), `include_items`, `n_item`, `factor_name`,
#'   `lev1`, `lev2`, `intercept`, `effect`, `family`, `resp_name`, and family parameters such
#'   as `sigma`, `shift`, `thresholds`, or `phi`; within-design random effects use
#'   `subj_int_sd`, `subj_slope_sd`, `subj_corr`, `item_int_sd`, `item_slope_sd`, and
#'   `item_corr`.
#' @return A design specification as a nested list, ready for [simulate_design()],
#'   [spec_json()], or the power and precision functions.
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
  if (p$family %in% c("gaussian", "shifted_lognormal")) { resp$sigma <- p$sigma; resp$round <- 4L }
  if (p$family == "shifted_lognormal") resp$shift <- p$shift
  if (p$family == "ordinal") resp$thresholds <- as.numeric(strsplit(gsub("\\s", "", p$thresholds), ",")[[1]])
  if (p$family == "beta") resp$phi <- p$phi
  spec$response <- resp
  spec
}

#' Serialise a design specification to pretty-printed JSON
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
spec_json <- function(spec) jsonlite::toJSON(spec, auto_unbox = TRUE, pretty = TRUE, digits = NA)

#' Generate a self-contained, reproducible R script from a specification
#'
#' Embed the specification as an R list literal (via `deparse`, which round-trips exactly), so
#' that the returned script reproduces the design without any external file. This turns a
#' design built in the no-code application into a reproducible script; the application's Verify
#' button runs that script in a clean R session and confirms that it reproduces the data
#' bit-for-bit.
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
    "spec <- ", paste(deparse(spec), collapse = "\n"), "\n\n",
    "data <- simulate_design(spec)              # analysis-ready data frame\n",
    "# write.csv(data, \"data.csv\", row.names = FALSE)\n",
    "# pow  <- power_mixed(spec, n_sims = 200)   # simulation-based power + Type S/M\n"
  )
}
