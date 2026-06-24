# Bridge to a Bayesian workflow. This translates a pilotr design spec into a ready-to-fit brms
# model, comprising the response family, the fixed and random-effects formula, and a
# weakly informative prior set. pilotr simulates the data, while brms (Stan) fits the
# confirmatory Bayesian model. The function emits code that can be copied into a script, and it
# requires neither brms nor Stan to be installed.

#' Derive a brms formula, family, and priors from a design spec
#'
#' @param spec a design spec (path or list).
#' @param prior_scale SD of the Normal prior on fixed main effects (standardized scale).
#' @param interaction_scale SD of the Normal prior on interaction terms (default prior_scale/2).
#' @return Invisibly, a list with elements `formula`, `family`, `priors`, and `code`; the
#'   `code` element (a ready-to-fit `brms` model) is also printed to the console.
#' @examples
#' spec <- build_spec(list(name = "d", seed = 1, design_kind = "within",
#'   include_items = TRUE, n_subject = 20, n_item = 12, factor_name = "cond",
#'   lev1 = "a", lev2 = "b", intercept = 6, effect = 0.05,
#'   subj_int_sd = 0.12, subj_slope_sd = 0.04, subj_corr = 0.2,
#'   item_int_sd = 0.08, item_slope_sd = 0.02, item_corr = -0.1,
#'   family = "shifted_lognormal", resp_name = "RT", sigma = 0.3, shift = 200))
#' bridge <- brms_bridge(spec)
#' bridge$formula
#' @export
brms_bridge <- function(spec, prior_scale = 0.5, interaction_scale = NULL) {
  if (is.character(spec)) spec <- load_spec(spec)
  if (is.null(interaction_scale)) interaction_scale <- prior_scale / 2

  family_map <- list(gaussian = "gaussian()", lognormal = "lognormal()",
                     shifted_lognormal = "shifted_lognormal()", bernoulli = "bernoulli()",
                     poisson = "poisson()", ordinal = "cumulative()")
  family <- family_map[[spec$response$family]]
  if (is.null(family)) stop("no brms family mapping for '", spec$response$family, "'")

  fixed_terms <- names(spec$fixed$coefficients)
  rs <- spec$random
  re_terms <- vapply(names(rs), function(g)
    sprintf("(%s | %s)", paste(c("1", names(rs[[g]]$slopes)), collapse = " + "), g),
    character(1))
  rhs <- paste(c(fixed_terms, re_terms), collapse = " + ")
  formula <- sprintf("%s ~ %s", spec$response$name, rhs)

  priors <- c('prior(normal(0, 2.5), class = "Intercept")')
  for (term in fixed_terms) {
    sc <- if (grepl(":", term, fixed = TRUE)) interaction_scale else prior_scale
    priors <- c(priors, sprintf('prior(normal(0, %s), class = "b", coef = "%s")', sc, term))
  }
  priors <- c(priors, 'prior(normal(0, 1), class = "sd")')        # half-normal (sd >= 0)
  has_cor <- any(vapply(rs, function(g) length(names(g$slopes)) > 0, logical(1)))
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
  cat(code, "\n")
  invisible(list(formula = formula, family = family, priors = priors, code = code))
}
