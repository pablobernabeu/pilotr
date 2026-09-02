# Regenerate the three cached results the power and precision vignettes ship.
#
# The mixed-model fits behind those vignettes cost far more than CRAN's check-time
# budget allows, so the vignettes read the results from disk instead of computing
# them. The chunks below are the ones the vignettes display, so whenever a
# replicate loop, a seed rule or a returned column changes, run this script from
# the package directory and commit the caches it writes.
#
# Install the package before running it, since with more than one worker the
# replicates are fitted in separate processes that load the installed namespace.
#
#   Rscript tools/regenerate-vignette-caches.R [workers]

args <- commandArgs(trailingOnly = TRUE)
workers <- if (length(args)) as.integer(args[[1]]) else 1L

library(pilotr)

# The specification of the power vignette. The precision vignette declares the
# same design with a smaller effect, hence the two builders.
priming_spec <- function(effect) {
  build_spec(list(
    name = "priming", seed = 1, design_kind = "within", include_items = TRUE,
    n_subject = 24, n_item = 18,
    factor_name = "condition", lev1 = "related", lev2 = "unrelated",
    intercept = 6, effect = effect,
    subj_int_sd = 0.12, subj_slope_sd = 0.04, subj_corr = 0.2,
    item_int_sd = 0.08, item_slope_sd = 0.02, item_corr = -0.1,
    family = "shifted_lognormal", resp_name = "RT", sigma = 0.3, shift = 200))
}

spec_power <- priming_spec(0.06)
spec_precision <- priming_spec(0.05)

# power-analysis.Rmd, 'Crossed mixed-effects designs'. The object is stored whole
# rather than as a table, so that the vignette prints what a live run prints.
pm <- power_mixed(spec_power, n_sims = 20, workers = workers)
saveRDS(pm, "vignettes/power-mixed-cache.rds", version = 2)

# power-analysis.Rmd, 'A power curve'.
curve <- power_curve_mixed(spec_power, subject_ns = c(8, 12, 16, 24, 32, 44, 56),
                           n_sims = 50, workers = workers)
write.csv(curve, "vignettes/power-curve-cache.csv", row.names = FALSE)

# precision-rope.Rmd, 'Sweeping sample size'.
prc <- precision_curve(spec_precision, focal = c(effect = 0.05),
                       subject_ns = c(15, 30, 60, 100, 140, 180, 220, 260),
                       rope = 0.02, n_sims = 200, workers = workers)
write.csv(prc, "vignettes/precision-curve-cache.csv", row.names = FALSE)
