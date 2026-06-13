# Headless test of the app's functional core (no Shiny UI). Confirms the GUI inputs build a
# valid portable spec, that the spec round-trips through JSON faithfully, and that simulate
# + power work on it -- i.e. the no-code front-end drives the same engine as R and Python.

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
source(file.path(here, "spec_builder.R"))
for (f in c("core.R", "simulate.R", "power.R"))
  source(file.path(here, "..", "r", "simdgp", "R", f))

ok <- TRUE
check <- function(cond, msg) { cat(if (cond) "  [PASS] " else "  [FAIL] ", msg, "\n", sep = ""); ok <<- ok && cond }

# ---- 1. Between-subjects Gaussian (default inputs) ----
p_between <- list(name = "ui_between", seed = 2024, n_subject = 64, include_items = FALSE,
                  n_item = 24, design_kind = "between", factor_name = "group",
                  lev1 = "control", lev2 = "treatment", intercept = 100, effect = 5,
                  family = "gaussian", resp_name = "", sigma = 10)
spec <- build_spec(p_between)
d <- simulate_design(spec)
cat("=== between-subjects gaussian ===\n")
check(nrow(d) == 64, "64 rows")
check(all(c("subject", "group", "score") %in% names(d)), "expected columns present")
check(is.null(spec$random), "no random effects for between design")

# round-trip through JSON (the artifact a user downloads)
tmp <- tempfile(fileext = ".json"); writeLines(spec_json(spec), tmp)
d2 <- simulate_design(load_spec(tmp))
check(isTRUE(all.equal(d$score, d2$score)), "JSON round-trip reproduces identical data")

res <- power_design(spec, n_sims = 500)
check(res$power > 0.3 && res$power < 0.7, sprintf("power in plausible range (%.3f)", res$power))

# ---- 2. Within / crossed mixed-effects RT ----
p_within <- list(name = "ui_crossed", seed = 90210, n_subject = 30, include_items = TRUE,
                 n_item = 24, design_kind = "within", factor_name = "condition",
                 lev1 = "related", lev2 = "unrelated", intercept = 6.0, effect = 0.05,
                 subj_int_sd = 0.12, subj_slope_sd = 0.04, subj_corr = 0.2,
                 item_int_sd = 0.08, item_slope_sd = 0.02, item_corr = -0.1,
                 family = "shifted_lognormal", resp_name = "", sigma = 0.30, shift = 200)
spec2 <- build_spec(p_within)
d3 <- simulate_design(spec2)
cat("=== within crossed mixed-effects RT ===\n")
check(nrow(d3) == 30 * 24 * 2, "1440 rows (subjects x items x condition)")
check(!is.null(spec2$random$subject$slopes), "by-subject random slope present")
check(!is.null(spec2$random$item$slopes), "by-item random slope present")
check(min(d3$RT) > 200, "RTs above the 200 ms shift")

cat("\nGenerated spec (within crossed):\n"); cat(spec_json(spec2), "\n")
cat(if (ok) "\nALL APP-LOGIC CHECKS PASSED\n" else "\nSOME CHECKS FAILED\n")
quit(status = if (ok) 0 else 1)
