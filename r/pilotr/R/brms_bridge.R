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
#
# brms is not a dependency, so nothing in the package's tests can ask brms whether it accepts
# what is emitted here. tools/brms/check_bridge.R does that for every shipped example, and it is
# the check to run after changing this file.

# The families whose coefficients are on the scale of the response itself, or of its logarithm
# for the two lognormal families. A prior of fixed width means something different on every
# such response, so these priors are scaled to the data. The link families' coefficients are on
# the logit or log scale, which does not depend on the response's units, so their priors keep a
# fixed width per unit of each column.
.bridge_scaled_families <- c("gaussian", "lognormal", "shifted_lognormal", "exgaussian")

# A prior width as written into the emitted code. Four significant digits is finer than any
# prior needs and short enough to read. sprintf() writes a full stop whatever the session's
# OutDec, so the code parses in every locale, and as.numeric() lets an integer scale such as 1L
# through a format meant for doubles.
.prior_width <- function(x) sprintf("%.4g", as.numeric(x))

# An interaction key's parts as a set, so that "cond:age" and "age:cond" compare equal.
.part_set <- function(keys) {
  vapply(strsplit(as.character(keys), ":", fixed = TRUE),
         function(p) paste(sort(p, method = "radix"), collapse = ":"), character(1))
}

# The names brms gives the fixed coefficients, keyed by the specification's names.
#
# brms builds its design matrix through stats::terms(), which names an interaction after its
# parts in the order they first appear in the formula. A specification that lists `age` before
# `cond` therefore has a coefficient brms calls "age:cond", whatever order the key "cond:age"
# gives, and a prior or a hypothesis written for "cond:age" names nothing brms has, so brm()
# stops. Each key maps to the term label with the same set of parts. A key that R cannot read as
# a formula term keeps its own spelling, since the formula that carries it cannot work either.
.brms_labels <- function(keys) {
  keys <- as.character(keys)
  if (!length(keys)) return(stats::setNames(character(0), character(0)))
  labs <- tryCatch(attr(stats::terms(stats::reformulate(keys)), "term.labels"),
                   error = function(e) keys)
  hit <- match(.part_set(keys), .part_set(labs))
  stats::setNames(ifelse(is.na(hit), keys, labs[hit]), keys)
}

# The standard deviation of each fixed coefficient's prior, and of the random-effect standard
# deviations' prior, before formatting.
#
# For the scaled families, a coefficient's prior is normal(0, s * sd_y / sd_x), with sd_y the
# response's total standard deviation on the model's scale and sd_x the standard deviation of the
# coefficient's column. `s` is then the prior standard deviation of the standardised coefficient,
# which is how rstanarm scales its default priors in a Gaussian model. The unit-scale prior this
# replaces held a five-point effect near zero on a response with a residual SD of 10, whatever the
# data said, and the Savage-Dickey Bayes factor came out near 1.
#
# sd_y comes from response_variance(), which measures the variance on the realised design and on
# the scale the coefficients act on, averaging over the random effects exactly. sd_x is measured
# on the realised column, interactions included, because the SD of a product is the product of
# the SDs only when the parts are centred and independent, and a predictor need not be either.
#
# The columns are read off model_data(spec, simulate_design(spec)), the data set the emitted call
# asks for, so each is the column brms is given: the observed value of a predictor measured with
# error, and the product of the observed parts for an interaction. One simulation serves every
# column. Taking each column from .design_column() instead costs one simulation per coefficient,
# about five seconds on reading_time_continuous, and returns a constant column with the 1e-12
# jitter of its near-zero residual in place of an SD of 0.
.bridge_prior_widths <- function(spec, keys, prior_scale, interaction_scale) {
  s <- ifelse(grepl(":", keys, fixed = TRUE), interaction_scale, prior_scale)
  if (!spec[["response"]][["family"]] %in% .bridge_scaled_families)
    return(list(b = s, sd = 1))
  sd_y <- sqrt(response_variance(spec)[["total"]])
  d <- model_data(spec, simulate_design(spec, validate = FALSE))
  sd_x <- vapply(keys, function(k) {
    x <- d[[.us(k)]]
    v <- stats::sd(x)
    # A spread below all.equal()'s tolerance on the column's own magnitude is rounding error,
    # not variation, and NA marks a design with a single row.
    if (is.finite(v) && v > sqrt(.Machine$double.eps) * max(1, abs(x))) v else NA_real_
  }, numeric(1))
  # A column that does not vary has no SD to divide by, and an infinite width is not code brms
  # can read. Its coefficient is confounded with the intercept, which the user should hear about.
  flat <- is.na(sd_x)
  for (k in keys[flat])
    warning(sprintf(
      "'%s' does not vary in this design, so its coefficient cannot be told apart from the intercept, and its prior is scaled as if the column had a standard deviation of 1",
      k), call. = FALSE)
  sd_x[flat] <- 1
  list(b = s * sd_y / sd_x, sd = sd_y)
}

#' Derive a brms formula, family, and priors from a design spec
#'
#' @details
#' The formula names the analysis columns of the design: each factor's numeric contrast columns
#' and the predictors, with every interaction left for brms to form from its parts. A data set
#' as [simulate_design()] returns it holds the factors' level labels and not their contrast
#' columns, so it goes through [model_data()] first, as in
#' `model_data(spec, simulate_design(spec))`. The emitted call says so beside its `data`
#' argument. Collected data laid out in the same way go through `model_data()` too.
#'
#' brms names an interaction coefficient after its parts in the order they first appear in the
#' formula, as [stats::terms()] does, so a key `cond:age` in a specification that lists `age`
#' first is the coefficient `age:cond`. The formula and the priors use the name brms uses. A
#' specification that gives one interaction two keys, such as `cond:age` and `age:cond`, is
#' refused, since brms estimates a single coefficient for it.
#'
#' For `gaussian` and `exgaussian` responses, the coefficients are on the response's own scale,
#' and for `lognormal` and `shifted_lognormal` on the scale of its logarithm, so their priors are
#' scaled to the data the design produces. Take `sd_y` to be the square root of the total
#' variance that [response_variance()] reports, and `sd_x` the standard deviation of a
#' coefficient's column in `model_data(spec, simulate_design(spec))`, an interaction's product
#' column included. The coefficient's prior is then `normal(0, s * sd_y / sd_x)`, where `s` is
#' `prior_scale` for a main effect and `interaction_scale` for an interaction, and the
#' random-effect standard deviations take a half-normal prior with standard deviation `sd_y`.
#' For `bernoulli`, `poisson`, `ordinal` and `beta`, the coefficients are on the logit or log
#' scale, which does not depend on the response's units. A coefficient's prior is then
#' `normal(0, s)` per unit of its column, and the random-effect standard deviations take
#' `normal(0, 1)`. These widths suit contrast columns and standardised predictors, while a
#' predictor on a wider scale, such as age in years, gets a correspondingly wider prior per
#' standard deviation. Widths are written to four significant digits.
#'
#' A prior on the random-effect standard deviations is emitted only when the design has random
#' effects, and an LKJ prior on their correlations only when some grouping factor's effects are
#' correlated. brms refuses a prior on a parameter the model does not contain. The intercept and
#' the family's own parameters, such as the residual standard deviation, take brms's default
#' priors. For the continuous families, brms centres the intercept's default prior on the median
#' of the response, or of its logarithm for the two lognormal families.
#'
#' A Bayes factor computed against these priors depends on their widths (Kass and Raftery,
#' 1995), so a different `prior_scale` gives a different Bayes factor from the same data. For a
#' continuous family at the default of 0.5, a standardised effect of 0.5 lies one prior standard
#' deviation from zero.
#'
#' @param spec a design spec (path or list).
#' @param prior_scale The standard deviation of the normal prior on a standardised main effect
#'   for the continuous families, that is on the coefficient multiplied by its column's standard
#'   deviation and divided by the response's. For `bernoulli`, `poisson`, `ordinal` and `beta`,
#'   it is the standard deviation of the prior on the coefficient itself, on the link scale.
#' @param interaction_scale The same for interaction terms. Defaults to half of `prior_scale`.
#' @references Kass, R. E. and Raftery, A. E. (1995). Bayes factors. \emph{Journal of the
#'   American Statistical Association}, 90(430), 773-795. \doi{10.1080/01621459.1995.10476572}
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
  # Validated, since the priors for the continuous families are measured on the design the
  # specification produces, and an invalid one cannot be simulated.
  spec <- .as_spec(spec)
  if (!.is_scalar_number(prior_scale) || prior_scale <= 0)
    stop("`prior_scale` must be a single positive number", call. = FALSE)
  if (is.null(interaction_scale)) interaction_scale <- prior_scale / 2
  if (!.is_scalar_number(interaction_scale) || interaction_scale <= 0)
    stop("`interaction_scale` must be a single positive number", call. = FALSE)

  family_map <- list(gaussian = "gaussian()", lognormal = "lognormal()",
                     shifted_lognormal = "shifted_lognormal()", bernoulli = "bernoulli()",
                     poisson = "poisson()", ordinal = "cumulative()",
                     # brms exgaussian(): mu is the mean, matching how the family is simulated.
                     exgaussian = "exgaussian()",
                     beta = "Beta()")  # brms Beta(): logit mu + precision phi, as simulated
  resp <- spec[["response"]]
  family <- family_map[[resp[["family"]]]]
  if (is.null(family)) stop("no brms family mapping for '", resp[["family"]], "'")

  keys <- as.character(names(spec[["fixed"]][["coefficients"]]))
  labels <- unname(.brms_labels(keys))
  # Two keys with the same parts, such as "cond:age" and "age:cond", pass validation and are
  # simulated as two terms on one product column, whose effects add. brms has one coefficient for
  # that interaction and refuses a second prior on it, so no model written from such a
  # specification would run.
  dup <- unique(labels[duplicated(labels)])
  if (length(dup))
    stop(paste(vapply(dup, function(l) {
      q <- sprintf("'%s'", keys[labels == l])
      sprintf("the keys %s and %s name one interaction, which brms estimates as the single coefficient '%s'",
              paste(q[-length(q)], collapse = ", "), q[length(q)], l)
    }, character(1)), collapse = "; "),
    "; write each interaction once, with the sum of its values", call. = FALSE)
  rs <- spec[["random"]]
  # `|` only for groups the specification actually correlates, `||` otherwise, matching what
  # simulate_design() generates. Emitting `|` unconditionally, together with an LKJ prior, told
  # Stan to estimate a correlation the process had fixed at zero.
  re_terms <- vapply(names(rs), function(g)
    sprintf("(%s %s %s)", paste(c("1", names(rs[[g]][["slopes"]])), collapse = " + "),
            .re_bar(rs[[g]]), g),
    character(1))
  # The fixed terms come first and in the specification's order, which is the order the labels
  # were taken in, so brms reads every interaction under the name the priors give it.
  rhs <- c(labels, re_terms)
  formula <- sprintf("%s ~ %s", resp[["name"]],
                     if (length(rhs)) paste(rhs, collapse = " + ") else "1")

  # No intercept prior: brms's default is centred on the response, where a fixed unit-scale
  # prior pulled the intercept of a response near 100 to about 4.
  widths <- .bridge_prior_widths(spec, keys, prior_scale, interaction_scale)
  priors <- sprintf('prior(normal(0, %s), class = "b", coef = "%s")',
                    .prior_width(widths$b), labels)
  # brms refuses a prior on a parameter the model does not contain, so the prior on the
  # random-effect SDs needs some random effect to apply to. It was emitted for every design,
  # and brm() stopped on each one without random effects. A normal prior on an SD is half-normal,
  # since brms bounds the SDs at zero.
  if (length(rs)) priors <- c(priors, sprintf('prior(normal(0, %s), class = "sd")',
                                              .prior_width(widths$sd)))
  # An LKJ prior likewise needs some group with a correlation matrix to put it on.
  has_cor <- any(vapply(rs, function(g)
    length(names(g[["slopes"]])) > 0 && .re_correlated(g), logical(1)))
  if (has_cor) priors <- c(priors, 'prior(lkj(2), class = "cor")')

  code <- paste0(
    "library(brms)\n",
    "fit <- brm(\n",
    "  ", formula, ",\n",
    "  data   = your_data,   # model_data(spec, simulate_design(spec)) adds the contrast columns\n",
    "  family = ", family, ",\n",
    if (length(priors)) paste0("  prior  = c(\n    ", paste(priors, collapse = ",\n    "),
                               "\n  ),\n"),
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
