# Power-vs-N curve for the crossed mixed-effects RT design (items fixed at 24).
# Writes build/power_curve_mixed.csv for the figure. Slow: ~ (#grid x n_sims) lmer fits.

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
source(file.path(here, "..", "R", "core.R"))
source(file.path(here, "..", "R", "simulate.R"))
source(file.path(here, "..", "R", "parallel.R"))
source(file.path(here, "..", "R", "power_mixed.R"))

BUILD <- file.path(here, "..", "..", "..", "build")
dir.create(BUILD, showWarnings = FALSE, recursive = TRUE)
spec <- load_spec(file.path(here, "..", "..", "..", "spec", "examples", "crossed_mixed_rt.json"))

grid <- c(10, 20, 30, 40, 50)
n_sims <- 60
cat(sprintf("Computing mixed-effects power curve over subjects = {%s}, %d sims each...\n",
            paste(grid, collapse = ", "), n_sims))
t0 <- Sys.time()
curve <- power_curve_mixed(spec, grid, n_sims = n_sims)
write.csv(curve, file.path(BUILD, "power_curve_mixed.csv"), row.names = FALSE)
cat(sprintf("done in %.0f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
print(curve)
