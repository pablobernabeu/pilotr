# Bridge to a Bayesian workflow. This translates a pilotr design spec into a ready-to-fit brms
# model, comprising the response family, the fixed and random-effects formula, and a
# weakly informative prior set. pilotr simulates the data, while brms (Stan) fits the
# confirmatory Bayesian model. The function returns code that can be copied into a script, and
# it requires neither brms nor Stan to be installed.
#
# The model is returned as a `pilotr_bridge` object rather than written to the console, so that
# calling the function costs nothing in output when the result is assigned or read by another
# function. The print method below is what puts the code on screen, so a bare call at the
# console still shows it. Printing from the function itself wrote to standard output even when
# nothing was there to read it, which CRAN's policy on console output asks a function not to do,
# and which forced the one internal caller to divert the output around itself.

#' Derive a brms formula, family, and priors from a design spec
#'
#' @param spec a design spec (path or list).
#' @param prior_scale SD of the Normal prior on fixed main effects (standardised scale).
#' @param interaction_scale SD of the Normal prior on interaction terms (default prior_scale/2).
#' @return An object of class `pilotr_bridge`: a list with elements `formula`, `family`,
#'   `priors`, and `code`, the last being a ready-to-fit `brms` model. The object is returned
#'   visibly, and [print.pilotr_bridge()] writes `code` to the console, so a bare call shows
#'   the model while an assignment stays silent.
#' @seealso [print.pilotr_bridge()] for the display, and [model_formula()] for the frequentist
#'   counterpart of the formula.
#' @examples
#' spec <- build_spec(list(name = "d", seed = 1, design_kind = "within",
#'   include_items = TRUE, n_subject = 20, n_item = 12, factor_name = "cond",
#'   lev1 = "a", lev2 = "b", intercept = 6, effect = 0.05,
#'   subj_int_sd = 0.12, subj_slope_sd = 0.04, subj_corr = 0.2,
#'   item_int_sd = 0.08, item_slope_sd = 0.02, item_corr = -0.1,
#'   family = "shifted_lognormal", resp_name = "RT", sigma = 0.3, shift = 200))
#' bridge <- brms_bridge(spec)          # silent
#' bridge$formula
#' bridge                               # prints the ready-to-fit model
#' @export
brms_bridge <- function(spec, prior_scale = 0.5, interaction_scale = NULL) {
  if (is.character(spec)) spec <- load_spec(spec)
  if (is.null(interaction_scale)) interaction_scale <- prior_scale / 2

  family_map <- list(gaussian = "gaussian()", lognormal = "lognormal()",
                     shifted_lognormal = "shifted_lognormal()", bernoulli = "bernoulli()",
                     poisson = "poisson()", ordinal = "cumulative()",
                     # brms exgaussian(): mu is the mean, matching how the family is simulated.
                     exgaussian = "exgaussian()",
                     beta = "Beta()")  # brms Beta(): logit mu + precision phi, as simulated
  family <- family_map[[spec$response$family]]
  if (is.null(family)) stop("no brms family mapping for '", spec$response$family, "'")

  fixed_terms <- names(spec$fixed$coefficients)
  rs <- spec$random
  # `|` only for groups the specification actually correlates, `||` otherwise, matching what
  # simulate_design() generates. Emitting `|` unconditionally, together with an LKJ prior, told
  # Stan to estimate a correlation the process had fixed at zero.
  re_terms <- vapply(names(rs), function(g)
    sprintf("(%s %s %s)", paste(c("1", names(rs[[g]]$slopes)), collapse = " + "),
            .re_bar(rs[[g]]), g),
    character(1))
  rhs <- paste(c(fixed_terms, re_terms), collapse = " + ")
  formula <- sprintf("%s ~ %s", spec$response$name, rhs)

  priors <- c('prior(normal(0, 2.5), class = "Intercept")')
  for (term in fixed_terms) {
    sc <- if (grepl(":", term, fixed = TRUE)) interaction_scale else prior_scale
    priors <- c(priors, sprintf('prior(normal(0, %s), class = "b", coef = "%s")', sc, term))
  }
  priors <- c(priors, 'prior(normal(0, 1), class = "sd")')        # half-normal (sd >= 0)
  # An LKJ prior only makes sense when some group has a correlation matrix to put it on. brms
  # rejects a prior on a parameter the model does not contain, so this has to track the bars.
  has_cor <- any(vapply(rs, function(g)
    length(names(g$slopes)) > 0 && .re_correlated(g), logical(1)))
  if (has_cor) priors <- c(priors, 'prior(lkj(2), class = "cor")')

  code <- paste0(
    "library(brms)\n",
    "fit <- brm(\n",
    "  ", formula, ",\n",
    "  data   = your_data,\n",
    "  family = ", family, ",\n",
    "  prior  = c(\n    ", paste(priors, collapse = ",\n    "), "\n  ),\n",
    "  chains = 4, iter = 4000, warmup = 2000, cores = 4,\n",
    "  control = list(adapt_delta = 0.95)\n)")
  out <- list(formula = formula, family = family, priors = priors, code = code)
  class(out) <- c("pilotr_bridge", "list")
  out
}

#' Print a brms bridge
#'
#' Writes the ready-to-fit `brms` model held in the object's `code` element, which is the form
#' the bridge is meant to be read in and the one to copy into a script. The other elements
#' (`formula`, `family`, `priors`) are the same model in parts, for a caller assembling its own
#' code, and are left to `str()` or to `$`.
#'
#' @param x A `pilotr_bridge` object, as returned by [brms_bridge()].
#' @param ... Ignored, present for consistency with the generic.
#' @return `x`, invisibly.
#' @examples
#' spec <- build_spec(list(name = "d", seed = 1, design_kind = "between",
#'   factor_name = "g", lev1 = "a", lev2 = "b", n_subject = 20,
#'   intercept = 0, effect = 0.4, family = "gaussian", resp_name = "y", sigma = 1))
#' print(brms_bridge(spec))
#' @export
print.pilotr_bridge <- function(x, ...) {
  # cat(), like print.pilotr_power(), so that the whole object reaches standard output as one
  # block: message() would put it on the stream knitr collects separately, splitting one
  # printed object across two boxes on the documentation site.
  cat(x$code, "\n")
  invisible(x)
}
