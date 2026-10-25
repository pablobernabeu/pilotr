# Portable numerical core for pilotr, mirrored operation for operation by python/pilotr/core.py.
# Shared L'Ecuyer (1988) combined LCG, Wichura AS 241 inverse-normal, hand-rolled
# Cholesky, and inverse-CDF response transforms. All LCG products stay below 2^53, so
# the arithmetic is exact in IEEE-754 doubles and matches the Python integer port. Values that
# pass through exp(), log() or pow() are bit-identical only within the limits spec/SPEC.md sets
# out under Scope of the guarantee.

.M1 <- 2147483563; .M2 <- 2147483399; .A1 <- 40014; .A2 <- 40692

#' Create a shared cross-language random-number generator
#'
#' Build the combined linear congruential generator (L'Ecuyer 1988) used by both the R and
#' Python implementations, so that a given specification and seed produce identical data in
#' either language wherever the two share a maths library. The uniforms are exact everywhere.
#' The draw-order contract, and where a value can differ in its last bits, are documented in the
#' specification at <https://github.com/pablobernabeu/pilotr/blob/main/spec/SPEC.md>.
#'
#' @param seed A single number used to seed the generator (coerced to a non-negative integer).
#' @return A list of three functions: `uniform()` returns one standard-uniform draw,
#'   `normal()` returns one standard-normal draw, and `normals(k)` returns a length-`k`
#'   numeric vector of standard-normal draws.
#' @examples
#' rng <- make_rng(1)
#' rng$uniform()
#' rng$normals(3)
#' @export
make_rng <- function(seed) {
  e <- new.env(parent = emptyenv())
  # Truncated, not rounded, to match `int(abs(seed))` in the Python twin. validate_spec refuses a
  # non-whole seed, so this only closes the validate = FALSE path the replicate loops use, where
  # rounding would have handed the two engines different data from one specification.
  seed <- abs(trunc(as.numeric(seed)))
  e$s1 <- 1 + (seed %% (.M1 - 1))
  e$s2 <- 1 + ((.A2 * e$s1) %% (.M2 - 1))
  unif <- function() {
    e$s1 <- (.A1 * e$s1) %% .M1
    e$s2 <- (.A2 * e$s2) %% .M2
    d <- e$s1 - e$s2
    if (d < 1) d <- d + (.M1 - 1)
    d / .M1
  }
  for (i in 1:10) unif()  # warm-up
  list(
    uniform = unif,
    normal  = function() as241(unif()),
    normals = function(k) vapply(seq_len(k), function(i) as241(unif()), numeric(1))
  )
}

#' Inverse normal cumulative distribution function (Wichura's AS 241)
#'
#' Compute the standard-normal quantile for a probability using Wichura's (1988) Algorithm
#' AS 241 (the PPND16 routine), the same algorithm underlying
#' [stats::qnorm()].
#'
#' The tail branch, which about 15% of uniform draws reach, calls `log()`. Where R and 'Python'
#' do not share a maths library, the two can therefore return different doubles for the same `p`.
#' On Windows, 14 of 1,000,000 draws differed, by 1 to 3 units in the last place.
#'
#' @param p A probability in the open interval (0, 1).
#' @return The standard-normal quantile (a numeric value) corresponding to `p`.
#' @examples
#' as241(0.975)   # about 1.959964
#' @export
as241 <- function(p) {
  q <- p - 0.5
  if (abs(q) <= 0.425) {
    r <- 0.180625 - q * q
    num <- (((((((2509.0809287301226727 * r + 33430.575583588128105) * r +
                 67265.770927008700853) * r + 45921.953931549871457) * r +
               13731.693765509461125) * r + 1971.5909503065514427) * r +
             133.14166789178437745) * r + 3.387132872796366608)
    den <- (((((((5226.495278852854561 * r + 28729.085735721942674) * r +
                 39307.89580009271061) * r + 21213.794301586595867) * r +
               5394.1960214247511077) * r + 687.1870074920579083) * r +
             42.313330701600911252) * r + 1.0)
    return(q * num / den)
  }
  r <- if (q < 0) p else 1.0 - p
  r <- sqrt(-log(r))
  if (r <= 5.0) {
    r <- r - 1.6
    num <- (((((((7.7454501427834140764e-4 * r + 0.0227238449892691845833) * r +
                 0.24178072517745061177) * r + 1.27045825245236838258) * r +
               3.64784832476320460504) * r + 5.7694972214606914055) * r +
             4.6303378461565452959) * r + 1.42343711074968357734)
    den <- (((((((1.05075007164441684324e-9 * r + 5.475938084995344946e-4) * r +
                 0.0151986665636164571966) * r + 0.14810397642748007459) * r +
               0.68976733498510000455) * r + 1.6763848301838038494) * r +
             2.05319162663775882187) * r + 1.0)
  } else {
    r <- r - 5.0
    num <- (((((((2.01033439929228813265e-7 * r + 2.71155556874348757815e-5) * r +
                 0.0012426609473880784386) * r + 0.026532189526576123093) * r +
               0.29656057182850489123) * r + 1.7848265399172913358) * r +
             5.4637849111641143699) * r + 6.6579046435011037772)
    den <- (((((((2.04426310338993978564e-15 * r + 1.4215117583164458887e-7) * r +
                 1.8463183175100546818e-5) * r + 7.868691311456132591e-4) * r +
               0.0148753612908506148525) * r + 0.13692988092273580531) * r +
             0.59983220655588793769) * r + 1.0)
  }
  z <- num / den
  if (q < 0) -z else z
}

# Inner product of the first `m` elements, written as an explicit IEEE-754 double fold.
#
# Neither language's built-in reduction is a plain double fold, and the two disagree with each
# other: base R's sum() accumulates in 80-bit long double on x86, while CPython's sum() has
# applied Neumaier compensation to float sequences since version 3.12. Both are more accurate
# than a naive fold, but they are more accurate in different ways, so an inner product of
# length three or more can land on different doubles in the two ports. Spelling the loop out
# keeps the arithmetic identical, which is what the cross-language guarantee needs.
.dot <- function(a, b, m) {
  s <- 0
  if (m > 0) for (k in seq_len(m)) s <- s + a[k] * b[k]
  s
}

# Lower Cholesky factor (Cholesky-Banachiewicz), mirroring core.py$cholesky.
#
# A negative pivot means the covariance implied by the specified standard deviations and
# correlations is not positive definite, so no random-effect distribution has those moments.
# Clamping the pivot at zero and carrying on would not fail, but it would not honour the
# specification either: the factor would generate random effects whose standard deviations are
# several times the requested ones. A correlation set that base R's chol() rejects is therefore
# an error, because the alternative is plausible-looking data from a process the user never
# described.
#
# A pivot of exactly zero is left alone. It is a useful case, and no kind of failure: it is
# what a slope with a standard deviation of zero produces, which is how a term is held fixed
# while the rest of the structure is kept intact. The tolerance below absorbs the rounding of
# an exactly-singular matrix, so a term that is negative only numerically, with no true
# negativity behind it, is clamped at zero.
.cholesky <- function(cov, label = NULL, cols = NULL) {
  n <- nrow(cov); L <- matrix(0, n, n)
  tol <- .Machine$double.eps * n * max(c(diag(cov), 1))
  for (i in 1:n) for (j in 1:i) {
    s <- .dot(L[i, ], L[j, ], j - 1)
    if (i == j) {
      d <- cov[i, i] - s
      if (d < -tol) .stop_not_pd(label, cols, i, d)
      L[i, j] <- sqrt(max(d, 0))
    } else {
      L[i, j] <- if (L[j, j] != 0) (cov[i, j] - s) / L[j, j] else 0
    }
  }
  L
}

# Reports which grouping factor failed, which random-effect column the factorisation reached,
# and how negative the pivot was. The offending column is more actionable than the smallest
# eigenvalue of the whole matrix, because it points at the correlations the user has to change,
# and it costs no extra numerics, so the R and Python ports can report it identically.
.stop_not_pd <- function(label, cols, i, d) {
  where <- if (!is.null(cols) && length(cols) >= i) sprintf(" at column '%s'", cols[i]) else
    sprintf(" at position %d", i)
  involved <- if (!is.null(cols)) sprintf(" Check the correlations among %s.",
                                          paste(sprintf("'%s'", cols[seq_len(i)]), collapse = ", ")) else ""
  stop(sprintf(
    "the random-effect covariance for '%s' is not positive definite: the Cholesky factorisation failed%s (pivot %g).%s No random-effect distribution has the requested standard deviations and correlations.",
    if (is.null(label)) "unknown group" else label, where, d, involved), call. = FALSE)
}

.matvec <- function(L, z) vapply(seq_len(nrow(L)), function(i) .dot(L[i, ], z, length(z)), numeric(1))

.inv_logit <- function(x) if (x >= 0) 1 / (1 + exp(-x)) else { e <- exp(x); e / (1 + e) }

# The largest mean the inverse-CDF walk serves. The walk starts from exp(-lam), which is
# subnormal once lam passes -log of the smallest normal double, 1022 * log(2). A subnormal
# keeps fewer significant bits, and from a mean of about 738 the cumulative distribution no
# longer sums to 1. Counts came out biased there (a mean of 744 drew a mean of 733.5), or the
# walk never reached the uniform, until exp(-lam) underflowed to zero near 745. The literal is
# the double nearest 1022 * log(2), and core.py compares against the same one, so both twins
# refuse the same means.
.POISSON_MAX_MEAN <- 708.3964185322641

# Both twins raise this text byte for byte. R prints an infinite mean as "Inf", which is what
# core.py's twin normalises its "inf" to.
.stop_poisson_mean <- function(lam) {
  stop(sprintf(
    "the poisson mean exp(eta) = %g is too large for the inverse-CDF sampler: beyond a mean of about 708.4, exp(-mean) is subnormal and the walk loses accuracy. Lower the poisson intercept or coefficients until the implied mean is at most 708.",
    lam), call. = FALSE)
}

# Each term is p_k = (p_{k-1} * lam) / k, multiplied before it is divided. spec/SPEC.md fixes
# that bracketing, because the other one lands on a different double for most terms. Within the
# cap, the cumulative sum reaches every uniform the generator draws long before the iteration
# cap, so the check after the loop only keeps an out-of-range u from being returned as a count.
.poisson_inv <- function(lam, u) {
  if (lam > .POISSON_MAX_MEAN) .stop_poisson_mean(lam)
  p <- exp(-lam); cum <- p; k <- 0
  while (u > cum && k < 1e6) { k <- k + 1; p <- p * lam / k; cum <- cum + p }
  if (u > cum) .stop_poisson_mean(lam)
  k
}

.ordinal_inv <- function(eta, thresholds, u) {
  for (k in seq_along(thresholds)) if (u <= .inv_logit(thresholds[k] - eta)) return(k)
  length(thresholds) + 1
}

# Marsaglia-Tsang Gamma(shape, scale=1) from the shared RNG, mirroring core.py's gamma_mt()
# operation for operation, with the one exception set out below for the boost.
#
# The cube is written as two multiplications rather than `^3`. R's `^` special-cases small
# integer exponents into repeated multiplication, whereas Python's `**` calls the library
# pow(); measured over 200,000 draws in this generator's range the two disagree on a third of
# inputs, by up to 6 ulp. Two explicit multiplies are exactly rounded in both languages, and
# they also decide the rejection step, so aligning them keeps the accept/reject sequence (and
# hence the number of draws consumed) identical across the two ports.
#
# The boost u^(1 / shape) is computed with `^`, which calls pow(), as Python's `**` does, unless
# 1 / shape is exactly 2, 3 or 4. From R 4.6.0, a cube or a fourth power is repeated
# multiplication, so a shape of exactly 1/3 or 1/4 can draw a different variate here than in
# core.py. spec/SPEC.md sets the case out under Scope of the guarantee.
.gamma_mt <- function(rng, shape) {
  if (shape < 1) return(.gamma_mt(rng, shape + 1) * rng$uniform()^(1 / shape))
  d <- shape - 1 / 3; cc <- 1 / sqrt(9 * d)
  repeat {
    x <- rng$normal(); t <- 1 + cc * x; v <- t * t * t
    if (v <= 0) next
    if (log(rng$uniform()) < 0.5 * x * x + d - d * v + d * log(v)) return(d * v)
  }
}

# The doubles nearest 1 and 0 inside the open interval, where a beta response of exactly 1 or 0
# is moved. 1 - 2^-53 is the largest double below 1, and 2^-1074 the smallest positive double, a
# subnormal. Both products are exact, and core.py forms the same two from float_info.
.BETA_HI <- 1 - .Machine$double.eps / 2
.BETA_LO <- .Machine$double.xmin * .Machine$double.eps

# Both twins raise these three texts byte for byte. The numbers go through %g in both, and
# core.py spells a value that is not finite as R does.
.stop_lognormal_overflow <- function(family, eta) {
  stop(sprintf(
    "the %s response exp(eta + sigma * z) overflows at a linear predictor of %g, so no response can be drawn. Lower the intercept or the coefficients.",
    family, as.numeric(eta)), call. = FALSE)
}

.stop_beta_shape <- function(eta, phi) {
  stop(sprintf(
    "the beta mean invlogit(eta) at a linear predictor of %g and precision phi = %g leaves a shape parameter of zero, so no proportion can be drawn. Keep eta between about -745 and 36.7 and phi large enough that mu * phi and (1 - mu) * phi stay positive.",
    as.numeric(eta), as.numeric(phi)), call. = FALSE)
}

.stop_beta_underflow <- function(eta, phi) {
  stop(sprintf(
    "both Gamma variates underflowed to zero at a linear predictor of %g and precision phi = %g, so the proportion is undefined. Raise phi.",
    as.numeric(eta), as.numeric(phi)), call. = FALSE)
}

# exp(w) for the two lognormal families, where w = eta + sigma * z. exp() overflows beyond about
# 709.78, and R returned Inf for every such row while Python's math.exp() raised a bare
# OverflowError. The refusal names eta, since exp(eta) itself can be finite in a row whose
# normal draw carries w past the limit.
.lognormal_exp <- function(w, family, eta) {
  val <- exp(w)
  if (!is.finite(val)) .stop_lognormal_overflow(family, eta)
  val
}

# A beta response for linear predictor eta and precision phi. It is x / (x + y), where x and y are
# Gamma variates with shapes mu * phi and (1 - mu) * phi, drawn in that order.
#
# A shape of zero leaves no distribution to draw from. invlogit(eta) rounds to 1 from an eta of
# about 36.74 and underflows to 0 below about -745.1, and a tiny phi can take either product to
# zero by itself. R's Gamma sampler then returned 0, so the response came out as exactly 0, 1 or
# NaN, while Python's raised ZeroDivisionError. When both shapes are tiny, both variates can also
# underflow to zero and leave 0 / 0. Both cases are refused before the division.
#
# Otherwise, the ratio rounds to exactly 1 whenever y is below about 2^-54 of x, as happens at
# realistic means near a ceiling, and it is exactly 0 when x underflows. brms's Beta() refuses
# both values, since they lie outside the family's support. Each becomes the nearest double inside
# the interval, 1 - 2^-53 or 2^-1074. Where both variates were representable, a ratio that rounded
# to 1 lay between 1 - 2^-53 and 1, so no double inside the interval is nearer the draw. The move
# is a comparison and draws nothing, so the stream is unchanged.
.beta_response <- function(rng, eta, phi) {
  mu <- .inv_logit(eta)
  a <- mu * phi; b <- (1 - mu) * phi
  if (a == 0 || b == 0) .stop_beta_shape(eta, phi)
  x <- .gamma_mt(rng, a)
  y <- .gamma_mt(rng, b)
  if (x + y == 0) .stop_beta_underflow(eta, phi)
  p <- x / (x + y)
  if (p == 1) .BETA_HI else if (p == 0) .BETA_LO else p
}

# A beta response rounded to `ndp` decimals by R's round(), with a result of exactly 0 or 1 moved
# inside (0, 1). Rounding sends every value below half a unit in the last place kept to 0. That
# happened to 17% of the rows of one realistic low-mean design, and brms's Beta() refuses such
# data. The bounds are the nearest values inside the interval at that many decimals, lo =
# 1 / 10^ndp and hi = 1 - lo: round() itself returns lo up to 22 places and hi up to 15. 10^ndp is
# built by multiplication, exact up to 22 places and free of either language's pow(). Beyond 22
# places lo falls back to 2^-1074, and from 17 places, where 1 - lo rounds to 1, hi falls back to
# 1 - 2^-53, the bounds of an unrounded response.
.beta_round <- function(val, ndp) {
  val <- round(val, ndp)
  if (val > 0 && val < 1) return(val)
  dig <- floor(ndp + 0.5)  # the number of places round() itself uses
  lo <- .BETA_LO
  if (dig >= 1 && dig <= 22) {
    p10 <- 1
    for (i in seq_len(dig)) p10 <- p10 * 10
    lo <- 1 / p10
  }
  if (val <= 0) return(lo)
  hi <- 1 - lo
  if (hi < 1) hi else .BETA_HI
}
