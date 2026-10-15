# Variance decomposition of the linear predictor, and calibration of the response scale.
#
# A design analysis usually needs the outcome on a known scale. A region of practical equivalence
# or a smallest effect size of interest is normally stated in standard-deviation units, and that
# is only meaningful once the total variance of the outcome is known, which for a crossed design
# with several random-effect terms is not something a user can read off the specification.
#
# The decomposition is computed from the realised design. An algebraic version would need the
# variance of every contrast column and predictor, and the
# covariance of every interaction, under assumptions about independence and about how between-unit
# factors divide the units. Those assumptions fail in exactly the awkward cases, such as a
# partially crossed design or an unbalanced between-unit split, so the components are measured on
# the design the specification actually produces.

# The per-row linear predictor, exactly as the generative engine computes it.
#
# eta does not depend on the response family: the family only maps eta to an outcome, and every
# draw eta depends on (predictors, then the random effects) is consumed before the per-row
# response draw. Substituting a Gaussian response with a negligible standard deviation therefore
# recovers eta itself for any family, including the discrete ones from which it could not
# otherwise be recovered, without reimplementing the row construction here.
.eta_draw <- function(spec) {
  s <- spec
  s[["response"]] <- list(family = "gaussian", name = ".eta", sigma = 1e-12)
  simulate_design(s, validate = FALSE)[[".eta"]]
}

# Set every random-effect standard deviation to zero while keeping every key in place.
#
# Keeping the keys matters: the number of normal deviates drawn is set by how many random-effect
# columns a grouping factor has, and their magnitudes do not enter, so zeroing the standard
# deviations leaves the RNG stream untouched and changes only the contribution to eta. Removing
# the keys would shorten the stream and shift everything downstream, which would make the
# components incomparable.
.zero_random <- function(spec, keep = character(0)) {
  for (g in names(spec[["random"]])) {
    if (g %in% keep) next
    spec[["random"]][[g]][["intercept_sd"]] <- 0
    for (k in names(spec[["random"]][[g]][["slopes"]]))
      spec[["random"]][[g]][["slopes"]][[k]] <- 0
  }
  spec
}

# The realised per-row values of one design column, whether a contrast column, a continuous
# predictor or an interaction between them.
#
# Recovered by making that column the only thing in the linear predictor: no intercept, a single
# coefficient of one, and every random-effect standard deviation zeroed. Neither coefficients nor
# standard deviations are drawn from the RNG, so this leaves the stream exactly where the real run
# has it, and the values returned are the ones that run will use, including the latent value of
# a predictor measured with error, in place of its observed value.
.design_column <- function(spec, key) {
  s <- .zero_random(spec)
  s[["fixed"]][["intercept"]] <- 0
  s[["fixed"]][["coefficients"]] <- stats::setNames(list(1), key)
  .eta_draw(s)
}

# How a family's residual variance behaves when the rest of the design is rescaled. This decides
# both what `calibrate_response()` can tune and how it has to solve for the scaling factor.
#
#   "free"     a residual standard deviation the specification sets, which scales with everything
#              else and can be solved for directly
#   "constant" a residual fixed by the link, which no rescaling can move
#   "eta"      a residual that depends on the linear predictor, so rescaling the design changes it
.residual_kind <- function(family) {
  switch(family,
    gaussian = , lognormal = , shifted_lognormal = , exgaussian = "free",
    bernoulli = , ordinal = "constant",
    poisson = , beta = "eta",
    "free")
}

.has_free_residual <- function(family) identical(.residual_kind(family), "free")

# The variance the response family adds on top of eta, on the scale eta lives on.
#
# Every family has one, so `total` is a complete decomposition for all eight families, including
# the four with no explicit `sigma`. What differs is the scale it is expressed on. For the
# lognormal families it is the log of the response, which is also the scale the auto-derived
# analysis model works on. For the link families it is the latent scale on which the linear
# predictor lives, which is the scale a coefficient, and therefore a region of practical
# equivalence, is stated on. This is the same "distribution-specific variance" used to compute the
# intraclass correlation and R-squared of a generalised mixed model (Nakagawa, Johnson and
# Schielzeth, 2017).
#
# Three of the four link cases are derived from pilotr's own generative process, and are exact
# for it:
#
#   bernoulli and ordinal  A row is drawn as 1[u < invlogit(eta)], equivalently logit(u) < eta, so
#                          the latent error is logit(u) with u uniform, which is a standard
#                          logistic variate of variance pi^2 / 3. The ordinal family compares the
#                          same uniform against cumulative thresholds, so it inherits it.
#   beta                   For Y ~ Beta(a, b), logit(Y) has variance trigamma(a) + trigamma(b)
#                          exactly. Confirmed against the package's own sampler across a range of
#                          shape pairs.
#
# The Poisson case is the one approximation. A log-scale variance cannot be measured directly,
# because a count of zero has no logarithm, so it is approximated. Of the three forms in the
# literature the trigamma one is used here, following Nakagawa et al.'s recommendation for the log
# link. The three agree closely once the mean exceeds about 5 and diverge sharply below it, where
# a log-scale variance is barely a meaningful quantity: at a mean of 0.5 the trigamma form gives
# 4.93 against 1.10 for the lognormal approximation. For a design dominated by rare counts, read
# the Poisson residual as an order of magnitude and no more.
#
# `eta` is the per-row linear predictor, needed only by the two families whose residual depends on
# the mean.
.residual_variance <- function(resp, eta = NULL) {
  sigma <- resp[["sigma"]]
  switch(resp[["family"]],
    gaussian          = sigma^2,
    lognormal         = sigma^2,
    shifted_lognormal = sigma^2,
    # Normal plus exponential, independent, so the variances add.
    exgaussian        = sigma^2 + resp[["beta"]]^2,
    # The latent logistic error of the threshold comparison.
    bernoulli         = pi^2 / 3,
    ordinal           = pi^2 / 3,
    # Averaged over the realised design, which keeps it consistent with how the fixed and
    # random components are computed. Evaluating at the mean would not.
    poisson           = mean(trigamma(exp(eta))),
    beta              = {
      phi <- .orelse(resp[["phi"]], 10)
      mu <- vapply(eta, .inv_logit, numeric(1))
      mean(trigamma(mu * phi) + trigamma((1 - mu) * phi))
    },
    NA_real_)
}

#' Variance components of the linear predictor
#'
#' Decompose the variance of a design's linear predictor into the part contributed by the fixed
#' effects, the part contributed by each grouping factor's random effects, and the residual
#' variance added by the response family. Useful for putting a region of practical equivalence or
#' a smallest effect size of interest on a known scale, and for seeing which term dominates a
#' design before committing to it.
#'
#' @details
#' The fixed component is the variance, across rows, of the linear predictor with every
#' random-effect standard deviation set to zero. It is read off the design the specification
#' actually produces, so it needs no assumption about how a between-unit factor divides the units
#' or about whether the predictors are independent.
#'
#' Each grouping factor's component is the average over rows of `x' Sigma x`, where `Sigma` is the
#' covariance the specification asks for and `x` collects that row's values of the intercept and
#' each random-slope column. This is exact for the realised design: it averages over the
#' random-effect distribution analytically, which matters because drawing from it and
#' estimating a variance from the drawn effects of, say, 30 subjects carries a sampling error of
#' around a quarter of the component itself, far too much to calibrate a region of practical
#' equivalence against.
#'
#' The cost is one simulation for the fixed part plus one per distinct random-slope column, which
#' for a typical design is a handful.
#'
#' Every component is on the scale the linear predictor lives on, which is what makes them
#' comparable. For `lognormal` and `shifted_lognormal` that is the log of the response, which is
#' also the scale the auto-derived analysis model works on. For `bernoulli`, `poisson`, `ordinal`
#' and `beta` it is the latent scale behind the link, so the residual is the distribution-specific
#' variance used to compute the intraclass correlation and R-squared of a generalised mixed model
#' (Nakagawa, Johnson and Schielzeth, 2017).
#'
#' Three of those four are derived from the process pilotr simulates and are exact for it. A
#' `bernoulli` row is drawn as `1[u < invlogit(eta)]`, equivalently `logit(u) < eta`,
#' so the latent error is a standard logistic variate of variance `pi^2 / 3`, and `ordinal` compares
#' the same uniform against cumulative thresholds and inherits it. For `beta`, `logit(Y)` has
#' variance `trigamma(a) + trigamma(b)` exactly.
#'
#' `poisson` is the one approximation, because a count of zero has no logarithm and a log-scale
#' variance cannot be measured directly. The trigamma form is used, following the recommendation for
#' the log link. It agrees closely with the alternatives once the mean exceeds about 5 and diverges
#' sharply below it, where a log-scale variance is barely meaningful: at a mean of 0.5 it gives 4.93
#' against 1.10 for the lognormal approximation. Read the Poisson residual as an order of magnitude
#' when counts are rare.
#'
#' The total includes the variance the fixed effects contribute. It is the denominator of an
#' R-squared, and of an effect or a region of practical equivalence stated in standard deviations
#' of the outcome, but it is not the denominator of Cohen's d. For designs in which participants
#' respond to stimuli, Westfall, Kenny and Judd (2014) define d as the difference between two
#' condition means divided by the standard deviation within a condition. The square of that
#' standard deviation sums the random intercept variances, the residual and the random slope
#' variances, each slope variance weighted by the squared contrast code. The fixed effects are left
#' out (see also Brysbaert and Stevens, 2018). In a balanced design with codes that sum to zero,
#' that is `total - fixed`. Covariates belong to the fixed part, so a d defined this way is
#' conditional on them.
#'
#' A published d therefore converts to a coefficient as `b = d * sqrt(total - fixed) / (c2 - c1)`,
#' where `c2 - c1` is the difference between the two contrast codes. For codes of -0.5 and 0.5,
#' that difference is 1. The d is on the scale of the linear predictor, so for the families with a
#' link it is a latent-scale d. [calibrate_response()] with `denominator = "random"` rescales a
#' design so that `total - fixed` is 1, after which the coefficient of a contrast coded -0.5 and
#' 0.5 is d.
#'
#' @param spec A design specification (path or list).
#' @return A named list of variance components: `fixed`, one entry per grouping factor,
#'   `residual`, and `total` (their sum).
#' @references
#' Brysbaert, M. and Stevens, M. (2018). Power analysis and effect size in mixed effects models: A
#'   tutorial. \emph{Journal of Cognition}, 1(1), 9. \doi{10.5334/joc.10}
#'
#' Nakagawa, S., Johnson, P. C. D. and Schielzeth, H. (2017). The coefficient of determination R2
#'   and intra-class correlation coefficient from generalized linear mixed-effects models revisited
#'   and expanded. \emph{Journal of the Royal Society Interface}, 14(134), 20170213.
#'   \doi{10.1098/rsif.2017.0213}
#'
#' Westfall, J., Kenny, D. A. and Judd, C. M. (2014). Statistical power and optimal design in
#'   experiments in which samples of participants respond to samples of stimuli. \emph{Journal of
#'   Experimental Psychology: General}, 143(5), 2020-2045. \doi{10.1037/xge0000014}
#' @examples
#' spec <- pilotr_example("crossed_mixed_rt")
#' response_variance(spec)
#'
#' # The coefficient that gives the related-unrelated contrast, coded -0.5 and 0.5, a published d
#' # of 0.2 on the log scale of this design.
#' v <- response_variance(spec)
#' 0.2 * sqrt(v$total - v$fixed) / (0.5 - (-0.5))
#'
#' # Every family reports a residual, including those whose outcome is discrete, so the components
#' # are a complete decomposition and their ratios read as the design's intraclass correlations.
#' v <- response_variance(pilotr_example("ordinal_likert_between"))
#' round(1 - v$residual / v$total, 3)   # share of latent variance that is structural
#' @seealso [calibrate_response()] to rescale a design to a target total variance.
#' @export
response_variance <- function(spec) {
  spec <- .as_spec(spec)
  rs <- spec[["random"]]
  groups <- names(rs)

  # stats::var() of a single value is NA, but a one-row design (validate_spec permits n = 1) has no
  # across-row variation at all, which is a different thing from an unknown amount of it.
  # Reporting 0 keeps `total` the sum of the parts it is documented to be, and leaves
  # calibrate_response() with a total it can solve against, where a quietly dropped term would
  # leave it with one it cannot.
  eta_fixed <- .eta_draw(.zero_random(spec))
  out <- list(fixed = if (length(eta_fixed) > 1L) stats::var(eta_fixed) else 0)

  # One simulation per distinct slope column, shared across the grouping factors that use it.
  slope_keys <- unique(unlist(lapply(groups, function(g) names(rs[[g]][["slopes"]]))))
  columns <- stats::setNames(lapply(slope_keys, function(k) .design_column(spec, k)), slope_keys)

  for (g in groups) {
    re <- .ranef(rs[[g]], g)
    Sigma <- re$L %*% t(re$L)
    # x collects the intercept (always 1) and this group's slope columns, per row.
    X <- cbind(1, if (length(re$cols) > 1L)
      do.call(cbind, columns[re$cols[-1L]]) else NULL)
    # mean over rows of x' Sigma x, which averages over the random-effect distribution exactly.
    M <- crossprod(X) / nrow(X)
    out[[g]] <- sum(Sigma * M)
  }

  # The full linear predictor, needed only by the families whose residual depends on the mean.
  resp <- spec[["response"]]
  eta <- if (identical(.residual_kind(resp[["family"]]), "eta")) .eta_draw(spec) else NULL
  out$residual <- .residual_variance(resp, eta)
  # A plain sum, with nothing silenced. validate_spec() admits only the families
  # .residual_variance() covers, so a missing component would be a fault worth surfacing rather
  # than one worth counting as zero and still calling the result the sum.
  out$total <- sum(vapply(out, identity, numeric(1)))
  out
}

# The variance calibrate_response() brings to the target, split into the structural part, which
# scales by k squared, and the residual. For denominator = "total", it is the whole total. For
# denominator = "random", it is the total less the fixed component, which is what Cohen's d
# divides by.
#
# The random part is summed from the grouping factors' components. Taking it as
# total - fixed - residual can leave a rounding remainder where a design with no random effects
# should give exactly zero. A positive remainder would then let a residual fixed by the link be
# closed by an enormous factor, where the design should be refused.
.calibrated_variance <- function(spec, parts, denominator) {
  if (identical(denominator, "total"))
    return(list(total = parts$total, structural = parts$total - parts$residual,
                residual = parts$residual, label = "fixed or random"))
  random <- sum(vapply(names(spec[["random"]]), function(g) parts[[g]], numeric(1)))
  list(total = random + parts$residual, structural = random, residual = parts$residual,
       label = "random-effect")
}

# The factor by which to multiply every scale-carrying term so that the calibrated variance `v`
# reaches the target, given how this family's residual responds to that multiplication.
#
# Scaling the coefficients and the random-effect standard deviations by k, and holding the
# intercept b0, moves each row's linear predictor from eta to b0 + k * (eta - b0). The design
# values are untouched, the random-effect covariance scales by k squared, and the number of draws
# consumed does not change, so the same deviates are drawn and the identity holds row by row.
# Scaling the intercept as well, as scale_intercept = TRUE does, is the same identity with b0 = 0,
# eta(k) = k * eta(1). Either is exact in real arithmetic. In floating point, a fresh simulation
# agrees with it to within 1e-9 of the linear predictor's standard deviation, on the shipped
# examples and for k from 0.1 to 10. That is far below anything a variance summary resolves. The
# identity is what lets the eta-dependent families be solved without simulating again at every
# candidate k.
#
# The intercept adds nothing to the fixed or random components, so the structural part always
# scales by k squared, and the three residual kinds differ only in what happens to the residual:
#
#   free      it scales too, so the target is reached in closed form
#   constant  it cannot move, so the structural part alone has to make up the difference, which is
#             impossible when the link's own residual already exceeds the target
#   eta       it moves with the linear predictor, so k is found numerically, cheaply, by reusing the
#             one eta vector
.calibration_factor <- function(spec, v, target_var, kind, scale_intercept = FALSE) {
  family <- spec[["response"]][["family"]]

  if (identical(kind, "free")) return(sqrt(target_var / v$total))

  if (identical(kind, "constant")) {
    if (target_var <= v$residual)
      stop(sprintf(
        "the %s family's own residual variance is %.6g on the latent scale, which already meets or exceeds the target of %.6g, and no rescaling of the design can move it; raise `target_var` above %.6g",
        family, v$residual, target_var, v$residual),
        call. = FALSE)
    if (v$structural <= 0)
      stop(sprintf("this design has no %s variance to rescale", v$label), call. = FALSE)
    return(sqrt((target_var - v$residual) / v$structural))
  }

  # kind == "eta": total(k) = k^2 * structural + residual(b0 + k * (eta - b0)), where b0 is the
  # intercept, or 0 when it is scaled with the rest. The residual is recomputed from the moved
  # linear predictor, with no fresh simulation. Both residuals are convex in eta (checked
  # numerically for eta from -12 to 12, and for the beta with phi from 0.2 to 500), so total(k) is
  # convex in k. At k = 0 every row sits at b0, and total(0) is the residual the intercept alone
  # implies. Centred contrasts and predictors make the rows' deviations from b0 average zero, apart
  # from the small mean that the realised predictor values and random effects leave. The total is
  # then flat at k = 0 and rises from there, so a target at or below total(0) is out of reach and
  # anything above it is bracketed and solved.
  # When the deviations do not average zero, the total can fall before it rises, and a target in
  # that dip is refused although some k reaches it. With the intercept held, that happens when a
  # contrast is coded 0 and 1 or a predictor's mean is not zero. Coded that way,
  # poisson_counts_between falls from 0.2213 at k = 0 to 0.2096 at k = 0.485. With b0 = 0, the
  # deviations average the intercept itself, and scale_intercept = TRUE keeps that refusal, as it
  # keeps the rest of the old behaviour. Both refusals say where the search starts.
  eta1 <- .eta_draw(spec)
  b0 <- if (scale_intercept) 0 else spec[["fixed"]][["intercept"]]
  total_at <- function(k)
    k * k * v$structural + .residual_variance(spec[["response"]], b0 + k * (eta1 - b0))
  lo <- total_at(0)
  if (target_var <= lo) {
    if (!scale_intercept)
      stop(sprintf(
        "with every coefficient and random-effect standard deviation set to zero, the intercept alone gives the %s family a residual variance of %.6g on the latent scale, which already meets or exceeds the target of %.6g; raise `target_var` above %.6g",
        family, lo, target_var, lo), call. = FALSE)
    # Holding the intercept starts the search from the residual at the intercept, so it is worth
    # suggesting only when that residual lies below the target.
    held <- .residual_variance(spec[["response"]], spec[["fixed"]][["intercept"]])
    stop(sprintf(
      "with every term set to zero, the intercept included, the %s family's residual variance at a linear predictor of zero is %.6g on the latent scale, which already meets or exceeds the target of %.6g; raise `target_var` above %.6g%s",
      family, lo, target_var, lo,
      if (held < target_var) ", or hold the intercept with scale_intercept = FALSE" else ""),
      call. = FALSE)
  }
  hi <- 1
  while (total_at(hi) < target_var && hi < 1e6) hi <- hi * 2
  if (total_at(hi) < target_var)
    stop("the target variance is not reachable by rescaling this design", call. = FALSE)
  stats::uniroot(function(k) total_at(k) - target_var, lower = 0, upper = hi,
                 tol = .Machine$double.eps^0.5)$root
}

#' Rescale a design to a target total variance
#'
#' Adjust a specification so that the total variance of its linear predictor, plus the residual
#' variance of its response family, comes to `target_var`. Calibrating to 1 puts the outcome on a
#' unit scale, which is what lets a region of practical equivalence or a smallest effect size of
#' interest be stated in standard-deviation units and read the same way across designs. The total
#' includes the variance the fixed effects contribute, so a coefficient calibrated against it is in
#' standard deviations of the outcome and is not Cohen's d. With `denominator = "random"`, the
#' target applies to the total less that fixed component, the variance that d divides by.
#'
#' @details
#' `tune = "sigma"` holds the fixed effects and the random-effect standard deviations where they
#' are and solves for the residual standard deviation. That keeps every effect size in the
#' specification as written, and is the right choice when those effects come from a pilot study or
#' from the literature. It fails when the structural variance already exceeds the target, since no
#' residual standard deviation can bring the total down, and the error says so rather than
#' returning a negative variance.
#'
#' `tune = "all"` holds the intercept and multiplies every coefficient and every random-effect
#' standard deviation, and the residual standard deviation where there is one, by a common factor.
#' For the families with a free residual that leaves every ratio between components unchanged, so
#' it rescales the outcome without altering the design's character. It is also the only option for
#' the families whose residual is fixed by the link.
#'
#' For those, the residual cannot be rescaled at all, so the structural part alone has to close the
#' gap and the target has to exceed the link's own contribution. A `bernoulli` or `ordinal` design
#' carries a latent residual of `pi^2 / 3`, about 3.29, so calibrating one to a total variance of 1
#' is not merely difficult but impossible, and the error says so rather than returning something
#' plausible. For `poisson` and `beta` the residual moves with the linear predictor, so the factor
#' is solved numerically; that costs one extra simulation rather than one per candidate, because
#' scaling the coefficients and the random-effect standard deviations by `k` moves each row's
#' linear predictor from `eta` to exactly `b0 + k * (eta - b0)`, where `b0` is the intercept. The
#' search starts at `k = 0`, where every row sits at the intercept, so a target at or below the
#' residual variance that the intercept alone implies is refused.
#'
#' The intercept is held because it is a location. It adds nothing to the fixed or random
#' components, so scaling it changes none of them and only moves the outcome. For `poisson` and
#' `beta`, whose residual depends on the location, scaling the intercept would reach the target
#' partly by moving the mean count or proportion. In the first example below, the median reaction
#' time of 600 ms stays near 600 ms. The intercept is the location when the contrasts and
#' predictors are centred, as codes of -0.5 and 0.5 are. Scaling the coefficient of a contrast
#' coded 0 and 1, or of a predictor whose mean is not zero, moves the location as well.
#'
#' `scale_intercept = TRUE` multiplies the intercept by the same factor, as pilotr 0.3.0 and 0.3.1
#' did. For `gaussian` and `exgaussian`, that is a change of units, and no power or precision
#' result for a coefficient changes. Through a log or logit link, it changes the outcome itself.
#' For the lognormal families, it raises the response, less any shift, to the power of the factor,
#' which turned that median of 600 ms into 16.6 hours. For `poisson`, it refuses targets that a
#' held intercept reaches, and it can return an intercept too large to simulate.
#'
#' `denominator` says what `target_var` is a target for. `"total"`, the default, is the total that
#' [response_variance()] reports, fixed component included. `"random"` is that total less the fixed
#' component, leaving the random-effect components and the residual. In a balanced design, that is
#' the variance within a condition. Westfall, Kenny and Judd (2014) define Cohen's d for designs in
#' which participants respond to stimuli by dividing a difference between condition means by the
#' square root of that variance. Once it is calibrated to 1, the coefficient of a two-level factor
#' coded -0.5 and 0.5 is d, on the scale of the linear predictor and conditional on any covariates.
#' The fixed component comes on top, so the total then exceeds `target_var`. With
#' `tune = "sigma"`, the residual standard deviation makes up what the random components leave of
#' the target. [response_variance()] gives the conversion of a published d into a coefficient.
#'
#' @param spec A design specification (path or list).
#' @param target_var The variance to calibrate to, as `denominator` defines it. Defaults to 1.
#' @param tune Either `"sigma"` (the default), which solves for the residual standard deviation
#'   alone, or `"all"`, which scales every variance-contributing term by a common factor.
#' @param scale_intercept With `tune = "all"`, whether the intercept is scaled with the rest.
#'   `FALSE` (the default) holds it, and `TRUE` reproduces pilotr 0.3.0 and 0.3.1.
#' @param denominator Either `"total"` (the default), the total variance with the fixed component
#'   included, or `"random"`, the random-effect components plus the residual, as in Cohen's d.
#' @return The specification, with the tuned parameters replaced.
#' @references
#' Brysbaert, M. and Stevens, M. (2018). Power analysis and effect size in mixed effects models: A
#'   tutorial. \emph{Journal of Cognition}, 1(1), 9. \doi{10.5334/joc.10}
#'
#' Westfall, J., Kenny, D. A. and Judd, C. M. (2014). Statistical power and optimal design in
#'   experiments in which samples of participants respond to samples of stimuli. \emph{Journal of
#'   Experimental Psychology: General}, 143(5), 2020-2045. \doi{10.1037/xge0000014}
#' @examples
#' spec <- pilotr_example("crossed_mixed_rt")
#' calibrated <- calibrate_response(spec, target_var = 1, tune = "all")
#' round(response_variance(calibrated)$total, 6)
#' calibrated$fixed$intercept   # held at 6, so the median reaction time stays near 600 ms
#'
#' # Calibrate the random components plus the residual, so that the coefficient of `cond`, coded
#' # -0.5 and 0.5, is Cohen's d on the log scale.
#' d_scale <- calibrate_response(spec, target_var = 1, tune = "all", denominator = "random")
#' d_scale$fixed$coefficients$cond
#' @seealso [response_variance()] for the decomposition this works from.
#' @export
calibrate_response <- function(spec, target_var = 1, tune = c("sigma", "all"),
                               scale_intercept = FALSE, denominator = c("total", "random")) {
  spec <- .as_spec(spec)
  tune <- match.arg(tune)
  denominator <- match.arg(denominator)
  if (!is.numeric(target_var) || length(target_var) != 1L || is.na(target_var) || target_var <= 0)
    stop("`target_var` must be a single positive number", call. = FALSE)
  if (!is.logical(scale_intercept) || length(scale_intercept) != 1L || is.na(scale_intercept))
    stop("`scale_intercept` must be TRUE or FALSE", call. = FALSE)

  parts <- response_variance(spec)
  v <- .calibrated_variance(spec, parts, denominator)
  family <- spec[["response"]][["family"]]
  beta <- spec[["response"]][["beta"]]
  kind <- .residual_kind(family)
  has_sigma <- identical(kind, "free")

  if (identical(tune, "sigma")) {
    if (!has_sigma)
      stop("the ", family, " family has no residual standard deviation to tune, ",
           "since its residual variance is fixed by the link rather than by a parameter; ",
           "use tune = \"all\" to rescale the rest of the design instead", call. = FALSE)
    random_only <- identical(denominator, "random")
    if (v$structural >= target_var)
      stop(sprintf(
        "the %s alone contribute a variance of %.6g, which already meets or exceeds the target of %.6g, so no residual standard deviation can reach it; lower the %s, raise `target_var`, or use tune = \"all\"",
        if (random_only) "random effects" else "fixed and random effects", v$structural,
        target_var, if (random_only) "random-effect standard deviations" else "effect sizes"),
        call. = FALSE)
    if (identical(family, "exgaussian")) {
      # The exponential component's variance is beta^2 and is being held fixed, so only the
      # normal component is free to absorb the difference.
      free <- target_var - v$structural - beta^2
      if (free <= 0)
        stop(sprintf(
          "with beta held at %.6g the exponential component alone contributes %.6g, leaving nothing for sigma; lower beta or use tune = \"all\"",
          beta, beta^2), call. = FALSE)
      spec[["response"]][["sigma"]] <- sqrt(free)
    } else {
      spec[["response"]][["sigma"]] <- sqrt(target_var - v$structural)
    }
    return(spec)
  }

  if (v$total <= 0)
    stop("this design has no variance to rescale", call. = FALSE)
  k <- .calibration_factor(spec, v, target_var, kind, scale_intercept)
  # The intercept is a location, which adds nothing to the fixed or random components, so it is
  # held unless the caller asks for the behaviour of pilotr 0.3.0 and 0.3.1.
  if (scale_intercept) spec[["fixed"]][["intercept"]] <- spec[["fixed"]][["intercept"]] * k
  for (nm in names(spec[["fixed"]][["coefficients"]]))
    spec[["fixed"]][["coefficients"]][[nm]] <- spec[["fixed"]][["coefficients"]][[nm]] * k
  for (g in names(spec[["random"]])) {
    re <- spec[["random"]][[g]]
    re[["intercept_sd"]] <- re[["intercept_sd"]] * k
    for (nm in names(re[["slopes"]])) re[["slopes"]][[nm]] <- re[["slopes"]][[nm]] * k
    spec[["random"]][[g]] <- re
  }
  if (has_sigma) spec[["response"]][["sigma"]] <- spec[["response"]][["sigma"]] * k
  if (identical(family, "exgaussian")) spec[["response"]][["beta"]] <- beta * k
  spec
}
