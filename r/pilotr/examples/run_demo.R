# Demo mirroring the Python run_demo.py: simulate both designs, recover parameters,
# run a simulation-based power + design analysis, and write CSVs for the parity check.

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
source(file.path(here, "..", "R", "core.R"))
source(file.path(here, "..", "R", "simulate.R"))
source(file.path(here, "..", "R", "parallel.R"))
source(file.path(here, "..", "R", "power.R"))

SPEC  <- file.path(here, "..", "..", "..", "spec", "examples")
BUILD <- file.path(here, "..", "..", "..", "build")
dir.create(BUILD, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Between-subjects Gaussian design ----
d <- simulate_design(file.path(SPEC, "between_2group_gaussian.json"))
write.csv(d, file.path(BUILD, "r_between.csv"), row.names = FALSE, quote = FALSE)
g0 <- d$score[d$group == "control"]; g1 <- d$score[d$group == "treatment"]
cat("=== between_2group_gaussian ===\n")
cat(sprintf("  N = %d ; control mean = %.3f, treatment mean = %.3f\n", nrow(d), mean(g0), mean(g1)))
cat(sprintf("  observed difference = %.3f (one noisy n=64 sample; true effect = 5.0)\n",
            mean(g1) - mean(g0)))

# ---- 2. Crossed mixed-effects reaction-time design ----
d2 <- simulate_design(file.path(SPEC, "crossed_mixed_rt.json"))
write.csv(d2, file.path(BUILD, "r_crossed.csv"), row.names = FALSE, quote = FALSE)
rel <- d2$RT[d2$condition == "related"]; unr <- d2$RT[d2$condition == "unrelated"]
cat("\n=== crossed_mixed_rt (subjects x items x condition) ===\n")
cat(sprintf("  N = %d rows ; mean RT related = %.1f ms, unrelated = %.1f ms\n",
            nrow(d2), mean(rel), mean(unr)))
cat(sprintf("  priming effect = %.1f ms (cond beta = 0.05 on log scale)\n", mean(unr) - mean(rel)))
cat("  first 3 rows:\n"); print(head(d2, 3))

# ---- continuous predictors + interactions + continuous random slopes ----
d3 <- simulate_design(file.path(SPEC, "reading_time_continuous.json"))
write.csv(d3, file.path(BUILD, "r_continuous.csv"), row.names = FALSE)
cat(sprintf("\n=== reading_time_continuous === N = %d rows; cols: %s\n", nrow(d3), paste(names(d3), collapse = ", ")))

# ---- additional grouping factor: subjects nested in higher-level clusters ----
d4 <- simulate_design(file.path(SPEC, "nested_clusters.json"))
write.csv(d4, file.path(BUILD, "r_nested.csv"), row.names = FALSE)
cat(sprintf("=== nested_clusters === N = %d rows; cols: %s\n", nrow(d4), paste(names(d4), collapse = ", ")))

# ---- Beta family (proportions) + partial crossing (item subset per subject) ----
for (nm in list(c("beta", "beta_proportion.json"), c("partial", "partial_crossing.json"))) {
  dd <- simulate_design(file.path(SPEC, nm[2]))
  write.csv(dd, file.path(BUILD, paste0("r_", nm[1], ".csv")), row.names = FALSE)
  cat(sprintf("=== %s === N = %d rows; cols: %s\n", nm[2], nrow(dd), paste(names(dd), collapse = ", ")))
}

# ---- 3. Simulation-based power + design analysis ----
res <- power_design(file.path(SPEC, "between_2group_gaussian.json"), n_sims = 2000)
cat("\n=== simulation-based power + design analysis ===\n")
cat(sprintf("  design: n=64 (32/group), d=0.5, alpha=.05, %d simulations\n", res$n_sims))
cat(sprintf("  power            = %.3f\n", res$power))
cat(sprintf("  Type S error     = %.4f\n", res$type_s))
cat(sprintf("  Type M (exagg.)  = %.3f\n", res$type_m))
cat(sprintf("  true effect = %.2f, mean estimate = %.3f\n", res$true_effect, res$mean_estimate))
