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
  s$response <- list(family = "gaussian", name = ".eta", sigma = 1e-12)
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
  for (g in names(spec$random)) {
    if (g %in% keep) next
    spec$random[[g]]$intercept_sd <- 0
    for (k in names(spec$random[[g]]$slopes)) spec$random[[g]]$slopes[[k]] <- 0
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
  s$fixed$intercept <- 0
  s$fixed$coefficients <- stats::setNames(list(1), key)
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
  switch(resp$family,
    gaussian          = resp$sigma^2,
    lognormal         = resp$sigma^2,
    shifted_lognormal = resp$sigma^2,
    # Normal plus exponential, independent, so the variances add.
    exgaussian        = resp$sigma^2 + resp$beta^2,
    # The latent logistic error of the threshold comparison.
    bernoulli         = pi^2 / 3,
    ordinal           = pi^2 / 3,
    # Averaged over the realised design, which keeps it consistent with how the fixed and
    # random components are computed. Evaluating at the mean would not.
    poisson           = mean(trigamma(exp(eta))),
    beta              = {
      phi <- if (is.null(resp$phi)) 10 else resp$phi
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
#' comparable and what makes their total the right denominator for an effect size. For
#' `lognormal` and `shifted_lognormal` that is the log of the response, which is also the scale the
#' auto-derived analysis model works on. For `bernoulli`, `poisson`, `ordinal` and `beta` it is the
#' latent scale behind the link, so the residual is the distribution-specific variance used to
#' compute the intraclass correlation and R-squared of a generalised mixed model (Nakagawa, Johnson
#' and Schielzeth, 2017).
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
#' @param spec A design specification (path or list).
#' @return A named list of variance components: `fixed`, one entry per grouping factor,
#'   `residual`, and `total` (their sum).
#' @references Nakagawa, S., Johnson, P. C. D. and Schielzeth, H. (2017). The coefficient of
#'   determination R2 and intra-class correlation coefficient from generalized linear mixed-effects
#'   models revisited and expanded. \emph{Journal of the Royal Society Interface}, 14(134),
#'   20170213. \doi{10.1098/rsif.2017.0213}
#' @examples
#' spec <- pilotr_example("crossed_mixed_rt")
#' response_variance(spec)
#'
#' # Every family reports a residual, including those whose outcome is discrete, so the components
#' # are a complete decomposition and their ratios read as the design's intraclass correlations.
#' v <- response_variance(pilotr_example("ordinal_likert_between"))
#' round(1 - v$residual / v$total, 3)   # share of latent variance that is structural
#' @seealso [calibrate_response()] to rescale a design to a target total variance.
#' @export
response_variance <- function(spec) {
  spec <- .as_spec(spec)
  groups <- names(spec$random)

  # stats::var() of a single value is NA, but a one-row design (validate_spec permits n = 1) has no
  # across-row variation at all, which is a different thing from an unknown amount of it.
  # Reporting 0 keeps `total` the sum of the parts it is documented to be, and leaves
  # calibrate_response() with a total it can solve against, where a quietly dropped term would
  # leave it with one it cannot.
  eta_fixed <- .eta_draw(.zero_random(spec))
  out <- list(fixed = if (length(eta_fixed) > 1L) stats::var(eta_fixed) else 0)

  # One simulation per distinct slope column, shared across the grouping factors that use it.
  slope_keys <- unique(unlist(lapply(groups, function(g) names(spec$random[[g]]$slopes))))
  columns <- stats::setNames(lapply(slope_keys, function(k) .design_column(spec, k)), slope_keys)

  for (g in groups) {
    re <- .ranef(spec$random[[g]], g)
    Sigma <- re$L %*% t(re$L)
    # x collects the intercept (always 1) and this group's slope columns, per row.
    X <- cbind(1, if (length(re$cols) > 1L)
      do.call(cbind, columns[re$cols[-1L]]) else NULL)
    # mean over rows of x' Sigma x, which averages over the random-effect distribution exactly.
    M <- crossprod(X) / nrow(X)
    out[[g]] <- sum(Sigma * M)
  }

  # The full linear predictor, needed only by the families whose residual depends on the mean.
  eta <- if (identical(.residual_kind(spec$response$family), "eta")) .eta_draw(spec) else NULL
  out$residual <- .residual_variance(spec$response, eta)
  # A plain sum, with nothing silenced. validate_spec() admits only the families
  # .residual_variance() covers, so a missing component would be a fault worth surfacing rather
  # than one worth counting as zero and still calling the result the sum.
  out$total <- sum(vapply(out, identity, numeric(1)))
  out
}

# The factor by which to multiply every scale-carrying term so that the total variance reaches the
# target, given how this family's residual responds to that multiplication.
#
# Scaling the intercept, the coefficients and the random-effect standard deviations all by k
# multiplies each row's linear predictor by k. The design values are untouched, the random-effect
# covariance scales by k squared, and the number of draws consumed does not change, so the same
# deviates are drawn and eta(k) is k times eta(1) row by row. That is exact in real arithmetic and
# holds to about 5e-13 in relative terms in floating point, which is far below anything a variance
# summary resolves. The identity is what lets the eta-dependent families be solved without
# simulating again at every candidate k.
#
# The structural part therefore always scales by k squared, and the three residual kinds differ only
# in what happens to the residual:
#
#   free      it scales too, so the target is reached in closed form
#   constant  it cannot move, so the structural part alone has to make up the difference, which is
#             impossible when the link's own residual already exceeds the target
#   eta       it moves with the linear predictor, so k is found numerically, cheaply, by reusing the
#             one eta vector
.calibration_factor <- function(spec, parts, target_var, kind) {
  structural <- parts$total - parts$residual

  if (identical(kind, "free")) return(sqrt(target_var / parts$total))

  if (identical(kind, "constant")) {
    if (target_var <= parts$residual)
      stop(sprintf(
        "the %s family's own residual variance is %.6g on the latent scale, which already meets or exceeds the target of %.6g, and no rescaling of the design can move it; raise `target_var` above %.6g",
        spec$response$family, parts$residual, target_var, parts$residual), call. = FALSE)
    if (structural <= 0)
      stop("this design has no fixed or random variance to rescale", call. = FALSE)
    return(sqrt((target_var - parts$residual) / structural))
  }

  # kind == "eta": total(k) = k^2 * structural + residual(k * eta), increasing in k, so bracket and
  # solve. The residual is recomputed from the scaled eta, with no fresh simulation.
  eta1 <- .eta_draw(spec)
  total_at <- function(k) k * k * structural + .residual_variance(spec$response, k * eta1)
  lo <- total_at(0)
  if (target_var <= lo)
    stop(sprintf(
      "even with every effect set to zero the %s family contributes a residual variance of %.6g on the latent scale, which already meets or exceeds the target of %.6g; raise `target_var` above %.6g",
      spec$response$family, lo, target_var, lo), call. = FALSE)
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
#' interest be stated in standard-deviation units and read the same way across designs.
#'
#' @details
#' `tune = "sigma"` holds the fixed effects and the random-effect standard deviations where they
#' are and solves for the residual standard deviation. That keeps every effect size in the
#' specification as written, and is the right choice when those effects come from a pilot study or
#' from the literature. It fails when the structural variance already exceeds the target, since no
#' residual standard deviation can bring the total down, and the error says so rather than
#' returning a negative variance.
#'
#' `tune = "all"` multiplies the intercept, every coefficient and every random-effect standard
#' deviation, and the residual standard deviation where there is one, by a common factor. For the
#' families with a free residual that leaves every ratio between components unchanged, so it
#' rescales the outcome without altering the design's character. It is also the only option for the
#' families whose residual is fixed by the link.
#'
#' For those, the residual cannot be rescaled at all, so the structural part alone has to close the
#' gap and the target has to exceed the link's own contribution. A `bernoulli` or `ordinal` design
#' carries a latent residual of `pi^2 / 3`, about 3.29, so calibrating one to a total variance of 1
#' is not merely difficult but impossible, and the error says so rather than returning something
#' plausible. For `poisson` and `beta` the residual moves with the linear predictor, so the factor
#' is solved numerically; that costs one extra simulation rather than one per candidate, because
#' scaling every term by `k` multiplies each row's linear predictor by exactly `k`.
#'
#' @param spec A design specification (path or list).
#' @param target_var The total variance to calibrate to. Defaults to 1.
#' @param tune Either `"sigma"` (the default), which solves for the residual standard deviation
#'   alone, or `"all"`, which scales every variance-contributing term by a common factor.
#' @return The specification, with the tuned parameters replaced.
#' @examples
#' spec <- pilotr_example("crossed_mixed_rt")
#' calibrated <- calibrate_response(spec, target_var = 1, tune = "all")
#' round(response_variance(calibrated)$total, 6)
#' @seealso [response_variance()] for the decomposition this works from.
#' @export
calibrate_response <- function(spec, target_var = 1, tune = c("sigma", "all")) {
  spec <- .as_spec(spec)
  tune <- match.arg(tune)
  if (!is.numeric(target_var) || length(target_var) != 1L || is.na(target_var) || target_var <= 0)
    stop("`target_var` must be a single positive number", call. = FALSE)

  parts <- response_variance(spec)
  kind <- .residual_kind(spec$response$family)
  has_sigma <- identical(kind, "free")

  if (identical(tune, "sigma")) {
    if (!has_sigma)
      stop("the ", spec$response$family, " family has no residual standard deviation to tune, ",
           "since its residual variance is fixed by the link rather than by a parameter; ",
           "use tune = \"all\" to rescale the rest of the design instead", call. = FALSE)
    structural <- parts$total - parts$residual
    if (structural >= target_var)
      stop(sprintf(
        "the fixed and random effects alone contribute a variance of %.6g, which already meets or exceeds the target of %.6g, so no residual standard deviation can reach it; lower the effect sizes, raise `target_var`, or use tune = \"all\"",
        structural, target_var), call. = FALSE)
    if (identical(spec$response$family, "exgaussian")) {
      # The exponential component's variance is beta^2 and is being held fixed, so only the
      # normal component is free to absorb the difference.
      free <- target_var - structural - spec$response$beta^2
      if (free <= 0)
        stop(sprintf(
          "with beta held at %.6g the exponential component alone contributes %.6g, leaving nothing for sigma; lower beta or use tune = \"all\"",
          spec$response$beta, spec$response$beta^2), call. = FALSE)
      spec$response$sigma <- sqrt(free)
    } else {
      spec$response$sigma <- sqrt(target_var - structural)
    }
    return(spec)
  }

  if (parts$total <= 0)
    stop("this design has no variance to rescale", call. = FALSE)
  k <- .calibration_factor(spec, parts, target_var, kind)
  spec$fixed$intercept <- spec$fixed$intercept * k
  for (nm in names(spec$fixed$coefficients))
    spec$fixed$coefficients[[nm]] <- spec$fixed$coefficients[[nm]] * k
  for (g in names(spec$random)) {
    spec$random[[g]]$intercept_sd <- spec$random[[g]]$intercept_sd * k
    for (nm in names(spec$random[[g]]$slopes))
      spec$random[[g]]$slopes[[nm]] <- spec$random[[g]]$slopes[[nm]] * k
  }
  if (has_sigma) spec$response$sigma <- spec$response$sigma * k
  if (identical(spec$response$family, "exgaussian"))
    spec$response$beta <- spec$response$beta * k
  spec
}
