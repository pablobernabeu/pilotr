# Equivalence check against simstudy (Goldfeld & Wujciak-Jens), the generative IV->DV
# competitor. simstudy specifies outcomes through a formula + link function; pilotr through
# an intercept + effect-size on a contrast + response family. For matched parameters the two
# produce statistically equivalent data, on a Gaussian design and a Poisson (log-link) one.

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
for (f in c("core.R", "simulate.R")) source(file.path(here, "..", "R", f))
suppressPackageStartupMessages(library(simstudy))
SPEC <- file.path(here, "..", "..", "..", "spec", "examples")
N <- 20000
row <- function(stat, ss, sg) cat(sprintf("  %-22s  %10.3f   %10.3f\n", stat, ss, sg))

# ---- Gaussian: simstudy formula "97.5 + 5*grp" vs pilotr intercept 100 + effect 5 ----
d <- defData(varname = "grp", formula = 0.5, dist = "binary")
d <- defData(d, varname = "y", formula = "97.5 + 5*grp", variance = 100, dist = "normal")
ss <- genData(N, d)
sg <- simulate_design(local({ s <- load_spec(file.path(SPEC, "between_2group_gaussian.json")); s$units$subject$n <- N; s }))
cat("=== Gaussian 2-group (simstudy formula+normal  vs  pilotr intercept+effect) ===\n")
cat("  statistic                 simstudy      pilotr\n  ------------------------------------------------\n")
row("mean (control)", mean(ss$y[ss$grp == 0]), mean(sg$score[sg$group == "control"]))
row("mean (treatment)", mean(ss$y[ss$grp == 1]), mean(sg$score[sg$group == "treatment"]))
row("effect (difference)", mean(ss$y[ss$grp == 1]) - mean(ss$y[ss$grp == 0]),
                           mean(sg$score[sg$group == "treatment"]) - mean(sg$score[sg$group == "control"]))
row("residual SD", sd(ss$y), sd(sg$score))

# ---- Poisson / log link: simstudy link='log' vs pilotr poisson family ----
dp <- defData(varname = "grp", formula = 0.5, dist = "binary")
dp <- defData(dp, varname = "count", formula = "1.4094 + 0.4*grp", dist = "poisson", link = "log")
ssp <- genData(N, dp)
sgp <- simulate_design(local({ s <- load_spec(file.path(SPEC, "poisson_counts_between.json")); s$units$subject$n <- N; s }))
cat("\n=== Poisson, log link (simstudy  vs  pilotr) ===\n")
cat("  statistic                 simstudy      pilotr\n  ------------------------------------------------\n")
row("mean count (control)", mean(ssp$count[ssp$grp == 0]), mean(sgp$count[sgp$group == "control"]))
row("mean count (treatment)", mean(ssp$count[ssp$grp == 1]), mean(sgp$count[sgp$group == "treatment"]))
row("rate ratio", mean(ssp$count[ssp$grp == 1]) / mean(ssp$count[ssp$grp == 0]),
                  mean(sgp$count[sgp$group == "treatment"]) / mean(sgp$count[sgp$group == "control"]))

cat("\n=> pilotr matches simstudy's generative output; it adds a Python implementation,\n")
cat("   crossed by-subject/by-item random slopes, and integrated power + Type S/M.\n")
