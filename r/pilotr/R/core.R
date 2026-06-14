# Portable numerical core for pilotr -- a bit-identical mirror of pilotr/core.py.
# Shared L'Ecuyer (1988) combined LCG, Wichura AS 241 inverse-normal, hand-rolled
# Cholesky, and inverse-CDF response transforms. All LCG products stay below 2^53, so
# the arithmetic is exact in IEEE-754 doubles and matches the Python integer port.

.M1 <- 2147483563; .M2 <- 2147483399; .A1 <- 40014; .A2 <- 40692

#' Create a shared cross-language RNG (see spec/SPEC.md for the contract).
#' @export
make_rng <- function(seed) {
  e <- new.env(parent = emptyenv())
  seed <- abs(round(as.numeric(seed)))
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

#' Wichura (1988) AS 241 inverse normal CDF (PPND16); same algorithm as R's qnorm.
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

# Lower Cholesky factor (Cholesky-Banachiewicz), mirroring core.py$cholesky.
.cholesky <- function(cov) {
  n <- nrow(cov); L <- matrix(0, n, n)
  for (i in 1:n) for (j in 1:i) {
    s <- if (j > 1) sum(L[i, 1:(j - 1)] * L[j, 1:(j - 1)]) else 0
    if (i == j) L[i, j] <- sqrt(max(cov[i, i] - s, 0))
    else        L[i, j] <- if (L[j, j] != 0) (cov[i, j] - s) / L[j, j] else 0
  }
  L
}

.matvec <- function(L, z) vapply(seq_len(nrow(L)), function(i) sum(L[i, ] * z), numeric(1))

.inv_logit <- function(x) if (x >= 0) 1 / (1 + exp(-x)) else { e <- exp(x); e / (1 + e) }

.poisson_inv <- function(lam, u) {
  p <- exp(-lam); cum <- p; k <- 0
  while (u > cum && k < 1e6) { k <- k + 1; p <- p * lam / k; cum <- cum + p }
  k
}

.ordinal_inv <- function(eta, thresholds, u) {
  for (k in seq_along(thresholds)) if (u <= .inv_logit(thresholds[k] - eta)) return(k)
  length(thresholds) + 1
}

# Marsaglia-Tsang Gamma(shape, scale=1) from the shared RNG; bit-identical with core.py.
.gamma_mt <- function(rng, shape) {
  if (shape < 1) return(.gamma_mt(rng, shape + 1) * rng$uniform()^(1 / shape))
  d <- shape - 1 / 3; cc <- 1 / sqrt(9 * d)
  repeat {
    x <- rng$normal(); v <- (1 + cc * x)^3
    if (v <= 0) next
    if (log(rng$uniform()) < 0.5 * x * x + d - d * v + d * log(v)) return(d * v)
  }
}
.beta_draw <- function(rng, a, b) { x <- .gamma_mt(rng, a); x / (x + .gamma_mt(rng, b)) }
