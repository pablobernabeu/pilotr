# Headless test of the app's functional core (no Shiny UI): GUI inputs -> portable spec ->
# JSON round-trip -> simulate + power. Confirms the no-code front-end drives the same engine
# as the R and Python packages.
#
# The engine is loaded exactly as the deployed browser build loads it: the subset of R/
# named in app-lite/engine-files.txt, and nothing else. Sourcing all of R/ here would keep
# this test green while a subset file grew a dependency on a file outside the subset, a
# break that then surfaces only in the deployed app.

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
manifest <- file.path(here, "..", "..", "..", "app-lite", "engine-files.txt")
engine <- trimws(sub("#.*$", "", readLines(manifest)))
engine <- engine[nzchar(engine)]
paths <- file.path(here, "..", "R", engine)
stopifnot(length(paths) > 0, all(file.exists(paths)))
for (f in paths) source(f)

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
# power_design() takes the Monte Carlo errors the app shows from mcse.R, so the subset has to
# carry that file as well.
check(all(c("power_mcse", "power_lo", "power_hi", "type_s_lo", "type_s_hi", "type_m_mcse")
          %in% names(res)) && res$power_lo < res$power && res$power < res$power_hi,
      "power_design() returns its Monte Carlo errors from the browser build's engine subset")
# The R script tab's script suggests the analysis that fits the design, by power_design()'s own
# rule in power.R, so the subset has to carry that file for spec_builder.R as well.
check(grepl("# pow <- power_design(spec, n_sims = 1000)", generate_r_script(spec), fixed = TRUE),
      "the exported script suggests power_design() from the browser build's engine subset")

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
check(grepl("# pow <- power_mixed(spec, n_sims = 200)", generate_r_script(spec2), fixed = TRUE),
      "the exported script of a crossed design suggests power_mixed()")

# Counterbalanced lists show each subject each item in one condition, rotated over subjects and
# items. They are a feature of specification version 0.4, which the engine subset reads.
spec3 <- build_spec(utils::modifyList(p_within, list(name = "ui_counterbalanced",
                                                     counterbalance = TRUE)))
d4 <- simulate_design(spec3)
cat("=== within crossed, counterbalanced lists ===\n")
check(identical(spec3$spec_version, "0.4") && isTRUE(spec3$factors[[1]]$counterbalance),
      "the counterbalanced option sets the field and declares 0.4")
check(nrow(d4) == 30 * 24 && anyDuplicated(d4[c("subject", "item")]) == 0,
      "720 rows (one per subject-item pair)")
tmp3 <- tempfile(fileext = ".json"); writeLines(spec_json(spec3), tmp3)
check(identical(simulate_design(load_spec(tmp3)), d4), "JSON round-trip reproduces identical data")

cat(if (ok) "\nALL APP-LOGIC CHECKS PASSED\n" else "\nSOME CHECKS FAILED\n")
quit(status = if (ok) 0 else 1)
