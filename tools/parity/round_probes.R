# Probe values for the rounding rule, each with the value R's own round() gives.
#
# R's round(x, digits) has followed one algorithm since R 4.0.0 (src/nmath/fround.c). It is the
# rounding rule of spec/SPEC.md, and the Python twin ports it as pilotr.core.fround(). CPython's
# built-in round() works on the exact binary value instead, and disagrees with R on a share of
# the doubles within a few units in the last place of a decimal tie. This script writes
# such doubles, random ones and the edge cases of every branch of fround.c, each with R's result,
# to python/tests/fixtures/r_round_probes.csv. python/tests/test_rounding.py checks the port
# against every row.
#
# Doubles are written in C's hexadecimal notation (sprintf("%a")), and every reader parses that
# exactly. A 17-digit decimal can be read one unit in the last place apart by a build of R whose
# long double is a plain double.
#
# Usage: Rscript tools/parity/round_probes.R [output-file]

if (getRversion() < "4.0.0")
  stop("R's rounding algorithm changed in 4.0.0, so run this script with R 4.0.0 or later",
       call. = FALSE)

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "..", ".."), mustWork = FALSE)
if (!nzchar(root) || is.na(root)) root <- normalizePath(".")
out <- if (length(args) >= 1) args[1] else
  file.path(root, "python", "tests", "fixtures", "r_round_probes.csv")

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(20261013)

# 1. Doubles within three units in the last place of a decimal tie: 1,250 at each of 1 to 4
#    digits, with integer parts from 1 to 10^7 and either sign.
near_tie <- function(d, n) {
  decade <- sample(0:6, n, replace = TRUE)
  k <- floor(runif(n, 10^(decade + d), 10^(decade + d + 1)))
  tie <- (k + 0.5) / c(10, 100, 1000, 10000)[d]   # the double nearest the decimal tie
  e <- floor(log2(tie))
  e <- e - (2^e > tie) + (2^(e + 1) <= tie)        # its binary exponent, whatever log2() rounds to
  x <- tie + sample(-3:3, n, replace = TRUE) * 2^(e - 52)
  data.frame(kind = "tie", x = sample(c(-1, 1), n, replace = TRUE) * x, digits = d)
}
ties <- do.call(rbind, lapply(1:4, near_tie, n = 1250))

# 2. Random doubles: 500 from a normal distribution with mean 100 and standard deviation 15 at 1 to
#    4 digits, and 500 spread over ten decades at 0 to 8 digits.
random <- rbind(
  data.frame(kind = "random", x = rnorm(500, 100, 15), digits = sample(1:4, 500, replace = TRUE)),
  data.frame(kind = "random", x = sample(c(-1, 1), 500, replace = TRUE) * 10^runif(500, -3, 7),
             digits = sample(0:8, 500, replace = TRUE)))

# 3. Ten doubles at each of 31 numbers of digits from 16 to 324. 10^d is exact only up to d = 22.
#    From d = 25, the binary powering R uses, R_pow_di(), and repeated multiplication can give
#    different powers, so a port has to power as R does. Beyond 308 digits, R splits the power in
#    two, and beyond 323 it returns x unchanged. Each x has 14 to 17 significant digits before the
#    rounding position, which takes in both sides of the shortcut that returns x unchanged once the
#    rounding position lies beyond its 15th significant digit.
many <- do.call(rbind, lapply(
  c(16:28, 30, 32, 35, 40, 60, 100, 200, 300, 305, 308:312, 315, 320, 323, 324),
  function(d) data.frame(kind = "digits",
                         x = runif(10, 1, 10) * 10^-(d - sample(13:16, 10, replace = TRUE)),
                         digits = d)))

# 4. The other branches of fround.c, with values on which R and CPython are known to differ and
#    values near either end of the double range. Digits of 0 give the nearest whole number, ties
#    to even. Digits that are not whole count as floor(digits + 0.5), without the shortcut for
#    exactly 0. Digits below -308 give 0.
edges <- list(
  c(0, 2), c(-0, 2), c(NaN, 1), c(Inf, 1), c(-Inf, 1), c(1.25, NaN), c(1.25, Inf), c(1.25, -Inf),
  c(0.5, 0), c(1.5, 0), c(2.5, 0), c(-0.5, 0), c(-2.5, 0), c(-0.3, 0), c(0.49999999999999994, 0),
  c(2^52 + 1, 0), c(2^50 + 0.5, 0), c(2^50 + 0.5, 0.25), c(1.23456, 2.5), c(1.5, -0.5),
  c(2.5, 0.375), c(1234.5678, -0.625),
  c(1234.5, -1), c(1250, -2), c(1350, -2), c(-1250, -2), c(123456789, -5), c(1e300, -300),
  c(1e300, -308), c(1e300, -309), c(-5, -400),
  c(1.05, 1), c(-1.05, 1), c(4.55, 1), c(-4.55, 1), c(394.255, 2), c(2.675, 2), c(1.005, 2),
  c(0.125, 2), c(0.375, 2), c(0.15, 1), c(0.25, 1), c(0.35, 1), c(0.45, 1), c(1.45, 1),
  c(123.456, 13), c(123.456, 14), c(0.1 + 0.2, 15), c(0.1 + 0.2, 16),
  c(5e-324, 2), c(5e-324, 323), c(5e-324, 324), c(1.5e-310, 315), c(-1.5e-310, 315),
  c(2.2250738585072014e-308, 310), c(.Machine$double.xmax, 2), c(-.Machine$double.xmax, -2))
edge <- data.frame(kind = "edge", x = vapply(edges, `[`, 0, 1), digits = vapply(edges, `[`, 0, 2))

probes <- rbind(ties, random, many, edge)
r <- round(probes$x, probes$digits)

lines <- c("kind,x,digits,r",
           paste(probes$kind, sprintf("%a", probes$x), sprintf("%.17g", probes$digits),
                 sprintf("%a", r), sep = ","))
con <- file(out, open = "wb")                     # binary: LF endings on every platform
writeLines(lines, con, sep = "\n")
close(con)
cat(sprintf("%d probes written to %s with %s\n", nrow(probes), out, R.version.string))
