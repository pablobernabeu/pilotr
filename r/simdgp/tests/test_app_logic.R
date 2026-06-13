# Headless test of the app's functional core (no Shiny UI): GUI inputs -> portable spec ->
# JSON round-trip -> simulate + power. Confirms the no-code front-end drives the same engine
# as the R and Python packages.

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
for (f in c("core.R", "simulate.R", "power.R", "spec_builder.R"))
  source(file.path(here, "..", "R", f))

ok <- TRUE
check <- function(cond, msg) { cat(if (cond) "  [PASS] " else "  [FAIL] ", msg, "\n", sep = ""); ok <<- ok && cond }

p_between <- list(name = "ui_between", seed = 2024, n_subject = 64, include_items = FALSE,
                  n_item = 24, design_kind = "between", factor_name = "group",
                  lev1 = "control", lev2 = "treatment", intercept = 100, effect = 5,
                  family = "gaussian", resp_name = "", sigma = 10)
spec <- build_spec(p_between); d <- simulate_design(spec)
cat("=== between-subjects gaussian ===\n")
check(nrow(d) == 64, "64 rows")
check(all(c("subject", "group", "score") %in% names(d)), "expected columns present")
check(is.null(spec$random), "no random effects for between design")
tmp <- tempfile(fileext = ".json"); writeLines(spec_json(spec), tmp)
check(isTRUE(all.equal(d$score, simulate_design(load_spec(tmp))$score)), "JSON round-trip reproduces identical data")
res <- power_design(spec, n_sims = 500)
check(res$power > 0.3 && res$power < 0.7, sprintf("power in plausible range (%.3f)", res$power))

p_within <- list(name = "ui_crossed", seed = 90210, n_subject = 30, include_items = TRUE,
                 n_item = 24, design_kind = "within", factor_name = "condition",
                 lev1 = "related", lev2 = "unrelated", intercept = 6.0, effect = 0.05,
                 subj_int_sd = 0.12, subj_slope_sd = 0.04, subj_corr = 0.2,
                 item_int_sd = 0.08, item_slope_sd = 0.02, item_corr = -0.1,
                 family = "shifted_lognormal", resp_name = "", sigma = 0.30, shift = 200)
spec2 <- build_spec(p_within); d3 <- simulate_design(spec2)
cat("=== within crossed mixed-effects RT ===\n")
check(nrow(d3) == 30 * 24 * 2, "1440 rows (subjects x items x condition)")
check(!is.null(spec2$random$subject$slopes) && !is.null(spec2$random$item$slopes), "by-subject AND by-item random slopes present")
check(min(d3$RT) > 200, "RTs above the 200 ms shift")

cat(if (ok) "\nALL APP-LOGIC CHECKS PASSED\n" else "\nSOME CHECKS FAILED\n")
quit(status = if (ok) 0 else 1)
