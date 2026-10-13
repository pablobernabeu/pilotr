# Monte Carlo uncertainty for the rates the design functions report.
#
# Every power and decision probability pilotr reports is a proportion over a finite number of
# replicates, so it carries its own uncertainty. Reporting it as a bare number invites the reader
# to treat it as exact. The defaults make this concrete: the curve functions default to 60
# replicates, at which a rate near 0.5 has a Monte Carlo standard error of 0.065, so a reported
# power of 0.80 is compatible with anything from roughly 0.70 to 0.90. That is not enough to
# support a claim about a design, and the only way a user can see it is if the number is shown.
#
# python/pilotr/mcse.py repeats every function here operation for operation, so that for the same
# counts the twins return the same doubles, and test-mcse.R and test_mcse.py pin them.

.Z95 <- 1.959963984540054   # qnorm(0.975), written out so the value cannot drift

# The plain Monte Carlo standard error of a proportion.
.mcse <- function(p, n) if (n <= 0) NA_real_ else sqrt(p * (1 - p) / n)

# Wilson score interval for a proportion.
#
# The plain standard error collapses to zero at p = 0 and p = 1, which reads as certainty
# precisely where there is least of it: nought significant results out of 60 is not evidence that
# power is exactly zero. The Wilson interval stays inside [0, 1] and keeps a sensible width at the
# boundaries, and it is the interval recommended over the Wald form for exactly this reason
# (Brown, Cai and DasGupta, 2001).
.wilson <- function(p, n, z = .Z95) {
  if (n <= 0) return(c(NA_real_, NA_real_))
  denom <- 1 + z * z / n
  centre <- (p + z * z / (2 * n)) / denom
  half <- z * sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
  lo <- centre - half; hi <- centre + half
  # At p = 0 and p = 1 the two terms cancel exactly in algebra but not in floating point, which
  # leaves a bound like 2.8e-17 where the answer is plainly 0. The endpoints are exact there, so
  # they are set by hand.
  if (p <= 0) lo <- 0
  if (p >= 1) hi <- 1
  c(max(0, lo), min(1, hi))
}

# A rate together with its Monte Carlo standard error and Wilson interval, as a flat named list
# suitable for splicing into a returned list or a data frame row.
.rate_with_error <- function(count, n, prefix) {
  p <- if (n > 0) count / n else NA_real_
  ci <- .wilson(if (is.na(p)) 0 else p, n)
  out <- list(p, .mcse(p, n), ci[1], ci[2])
  names(out) <- c(prefix, paste0(prefix, "_mcse"), paste0(prefix, "_lo"), paste0(prefix, "_hi"))
  out
}

# The Monte Carlo standard error of Type M: the standard deviation of the exaggeration ratios over
# the significant replicates, divided by the square root of their number. Type M averages over the
# significant replicates alone, often a handful, so it is the least settled number a design
# analysis reports. A standard deviation needs two values, so fewer give NA. The sums are plain
# double folds, written out because sum(), mean() and sd() accumulate in long double where the
# platform has one. The Python twin could not repeat those sums bit for bit, and nor could a
# platform without long double.
.type_m_mcse <- function(ratios) {
  n <- length(ratios)
  if (n < 2L) return(NA_real_)
  s <- 0
  for (r in ratios) s <- s + r
  m <- s / n
  ss <- 0
  for (r in ratios) ss <- ss + (r - m) * (r - m)
  sqrt(ss / (n - 1) / n)
}

# The Wilson interval for Type S, the share of the significant replicates whose estimate has the
# wrong sign. Its plain standard error is zero whenever no sign is wrong. That is the usual result,
# and the standard error would then claim certainty from what may be a handful of replicates.
# Type S is defined only against a non-zero true effect and over at least one significant
# replicate, so otherwise both bounds are NA.
.type_s_interval <- function(n_wrong, n_sig, true_effect) {
  if (is.na(true_effect) || true_effect == 0 || n_sig <= 0) return(c(NA_real_, NA_real_))
  .wilson(n_wrong / n_sig, n_sig)
}
