# Auto-derive the analysis model (data and lmer formula) from a design spec, so that design
# analysis and power can run from a spec alone, without a hand-coded formula. Interactions
# "a:b" become product columns named "a_b". Categorical factors become their numeric contrast
# columns. The response is log-transformed for the (shifted_)lognormal families. This targets
# the lmer-fittable families (gaussian, lognormal, shifted_lognormal). Other families require
# glmer with the appropriate link.

.us <- function(k) gsub(":", "_", k, fixed = TRUE)   # "a:b" -> "a_b"

#' Build the modelling data frame from a simulated data set + its spec.
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

#' Derive the lmer formula implied by a spec (response `.y` from `model_data`).
#' @export
model_formula <- function(spec) {
  fixed <- vapply(names(spec$fixed$coefficients), .us, character(1))
  re <- vapply(names(spec$random), function(g)
    sprintf("(%s | %s)", paste(c("1", vapply(names(spec$random[[g]]$slopes), .us, character(1))),
                               collapse = " + "), g), character(1))
  stats::as.formula(paste(".y ~", paste(c(fixed, re), collapse = " + ")))
}
