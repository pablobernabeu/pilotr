# Auto-derive the analysis model (data and lmer formula) from a design spec, so that design
# analysis and power can run from a spec alone, without a hand-coded formula. Interactions
# "a:b" become product columns named "a_b". Categorical factors become their numeric contrast
# columns. The response is log-transformed for the (shifted_)lognormal families. This targets
# the lmer-fittable families (gaussian, lognormal, shifted_lognormal). Other families require
# glmer with the appropriate link.

.us <- function(k) gsub(":", "_", k, fixed = TRUE)   # "a:b" -> "a_b"

#' Build the modelling data frame from a simulated data set and its specification
#'
#' Add the analysis response column `.y` (log-transformed for the lognormal families), the
#' numeric contrast columns implied by the categorical factors, and any interaction product
#' columns (an `a:b` coefficient becomes a column `a_b`).
#'
#' @param spec A design specification (path or list).
#' @param d A simulated data set, as returned by [simulate_design()].
#' @return The data frame `d` augmented with the `.y` response column and the contrast and
#'   interaction columns required by the auto-derived model.
#' @examples
#' spec <- build_spec(list(name = "d", seed = 1, design_kind = "between",
#'   factor_name = "group", lev1 = "a", lev2 = "b", n_subject = 20,
#'   intercept = 0, effect = 0.5, family = "gaussian", resp_name = "", sigma = 1))
#' head(model_data(spec, simulate_design(spec)))
#' @export
model_data <- function(spec, d) {
  resp <- spec$response; shift <- if (is.null(resp$shift)) 0 else resp$shift
  d$.y <- if (resp$family %in% c("lognormal", "shifted_lognormal")) log(d[[resp$name]] - shift) else d[[resp$name]]
  for (f in spec$factors) for (col in names(f$contrasts))      # contrast columns from labels
    d[[col]] <- f$contrasts[[col]][match(d[[f$name]], f$levels)]
  for (key in names(spec$fixed$coefficients)) if (grepl(":", key, fixed = TRUE)) {  # interaction products
    parts <- strsplit(key, ":", fixed = TRUE)[[1]]
    d[[.us(key)]] <- Reduce(`*`, lapply(parts, function(p) d[[p]]))
  }
  d
}

#' Derive the lmer formula implied by a specification
#'
#' Construct the mixed-model formula (with response `.y`, as produced by
#' [model_data()]) from the fixed-effect coefficients and random-effects
#' structure of a specification. Interaction coefficients written `a:b` become
#' formula terms `a_b`.
#'
#' @param spec A design specification (path or list).
#' @return A [stats::formula] object suitable for fitting with `lme4` or
#'   `lmerTest`.
#' @examples
#' spec <- build_spec(list(name = "d", seed = 1, design_kind = "within",
#'   include_items = TRUE, n_subject = 10, n_item = 8, factor_name = "cond",
#'   lev1 = "a", lev2 = "b", intercept = 6, effect = 0.05,
#'   subj_int_sd = 0.12, subj_slope_sd = 0.04, subj_corr = 0.2,
#'   item_int_sd = 0.08, item_slope_sd = 0.02, item_corr = -0.1,
#'   family = "shifted_lognormal", resp_name = "", sigma = 0.3, shift = 200))
#' model_formula(spec)
#' @export
model_formula <- function(spec) {
  fixed <- vapply(names(spec$fixed$coefficients), .us, character(1))
  re <- vapply(names(spec$random), function(g)
    sprintf("(%s | %s)", paste(c("1", vapply(names(spec$random[[g]]$slopes), .us, character(1))),
                               collapse = " + "), g), character(1))
  # Anchor the formula in the global environment rather than this function's evaluation
  # frame. print.formula omits its `<environment: ...>` line only for the global
  # environment, so any other choice trails a raw memory address that changes on every
  # run, which is noise in the printed result and makes the documented output impossible
  # to reproduce. Nothing is lost by it: every term is a column name that model_data()
  # creates, and those are always resolved from the `data` argument of the fit, never
  # through the formula's environment.
  stats::as.formula(paste(".y ~", paste(c(fixed, re), collapse = " + ")), env = globalenv())
}
