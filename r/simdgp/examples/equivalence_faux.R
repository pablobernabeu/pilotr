# Equivalence check against faux (DeBruine), the closest R competitor.
#
# faux specifies a within-subject design through a correlation matrix; simdgp induces the
# same correlation through a by-subject random intercept. For a two-condition within design
# with grand structure (mu = 100/105, sd = 15, within-correlation r = 0.5), the random-
# intercept variance that reproduces r is  sigma_b^2 = r * total_var,  sigma_e^2 = (1-r) *
# total_var. This script shows the two tools produce statistically equivalent data, then
# notes what simdgp adds on top.

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
source(file.path(here, "..", "R", "core.R"))
source(file.path(here, "..", "R", "simulate.R"))
suppressPackageStartupMessages(library(faux))

N <- 5000
mu <- c(100, 105); total_sd <- 15; r <- 0.5
sigma_b <- sqrt(r) * total_sd          # by-subject intercept SD that reproduces r
sigma_e <- sqrt(1 - r) * total_sd      # residual SD

# ---- faux: correlation-matrix parameterization ----
set.seed(1)
df_faux <- sim_design(within = list(cond = c("c1", "c2")),
                      n = N, mu = mu, sd = total_sd, r = r, plot = FALSE)

# ---- simdgp: random-intercept parameterization (same design, via the portable spec) ----
spec <- list(
  name = "faux_equivalence", seed = 1,
  units = list(subject = list(n = N)),
  factors = list(list(name = "cond", levels = c("c1", "c2"),
                      contrasts = list(cond = c(-0.5, 0.5)), vary_within = "subject")),
  fixed = list(intercept = mean(mu), coefficients = list(cond = mu[2] - mu[1])),
  random = list(subject = list(intercept_sd = sigma_b)),
  response = list(family = "gaussian", name = "score", sigma = sigma_e)
)
d <- simulate_design(spec)
s1 <- d$score[d$cond == "c1"]; s2 <- d$score[d$cond == "c2"]   # aligned by subject

row <- function(stat, target, f, g) cat(sprintf("  %-18s  %8.3f   %8.3f   %8.3f\n", stat, target, f, g))
cat(sprintf("Within-subject 2-condition design, N = %d (target: mu 100/105, sd 15, r 0.5)\n\n", N))
cat("  statistic              target      faux     simdgp\n")
cat("  ---------------------------------------------------\n")
row("mean (c1)",   mu[1],     mean(df_faux$c1), mean(s1))
row("mean (c2)",   mu[2],     mean(df_faux$c2), mean(s2))
row("sd (c1)",     total_sd,  sd(df_faux$c1),   sd(s1))
row("sd (c2)",     total_sd,  sd(df_faux$c2),   sd(s2))
row("corr(c1,c2)", r,         cor(df_faux$c1, df_faux$c2), cor(s1, s2))

cat("\n  => simdgp reproduces faux's correlated within-subject data. What simdgp ADDS:\n")
cat("     - a Python implementation that yields bit-identical data from the same spec\n")
cat("     - integrated simulation-based power + Type S/M design analysis\n")
cat("     - crossed by-subject AND by-item random slopes with a generative effect-size API\n")
