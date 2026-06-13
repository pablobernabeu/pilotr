# GUI-input -> portable-spec logic, as package functions (so the app stays a thin client
# and this logic is unit-testable). Pure R; no Shiny.

#' Default response-column name for a family.
#' @export
default_response_name <- function(family) {
  switch(family,
         gaussian = "score", shifted_lognormal = "RT",
         bernoulli = "accuracy", poisson = "count", ordinal = "rating", "outcome")
}

#' Build a simdgp design spec (a plain list) from a flat list of UI inputs.
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
  spec$response <- resp
  spec
}

#' Serialise a spec to pretty JSON (the portable artifact).
#' @export
spec_json <- function(spec) jsonlite::toJSON(spec, auto_unbox = TRUE, pretty = TRUE, digits = NA)
