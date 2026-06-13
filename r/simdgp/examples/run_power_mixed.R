# Simulation-based power + design analysis for the crossed mixed-effects RT design.

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
source(file.path(here, "..", "R", "core.R"))
source(file.path(here, "..", "R", "simulate.R"))
source(file.path(here, "..", "R", "power_mixed.R"))

spec <- file.path(here, "..", "..", "..", "spec", "examples", "crossed_mixed_rt.json")
n_sims <- 100

cat(sprintf("Fitting %d maximal mixed models (RT ~ cond + (1+cond|subject) + (1+cond|item))...\n", n_sims))
t0 <- Sys.time()
res <- power_mixed(spec, n_sims = n_sims)
cat(sprintf("done in %.0f s\n\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

cat("=== simulation-based power: crossed mixed-effects RT design ===\n")
cat(sprintf("  30 subjects x 24 items x 2 conditions; fixed effect (log scale) = %.3f\n", res$true_effect))
cat(sprintf("  converged fits   = %d / %d\n", res$n_converged, res$n_sims))
cat(sprintf("  power            = %.3f\n", res$power))
cat(sprintf("  Type S error     = %.4f\n", res$type_s))
cat(sprintf("  Type M (exagg.)  = %.3f\n", res$type_m))
cat(sprintf("  mean estimate    = %.4f (true %.3f)\n", res$mean_estimate, res$true_effect))
