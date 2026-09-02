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

# The R script tab and its download hand the reader a script that carries the design as a
# literal. Sourcing a file proves that it parses, not that its dependencies are staged, so the
# literal is evaluated here and simulated: that is what the app promises the script does.
cat("=== the R script the app exports ===\n")
exprs <- parse(text = generate_r_script(spec))
assign_spec <- Filter(function(e) is.call(e) && identical(e[[1]], as.name("<-")) &&
                        identical(e[[2]], as.name("spec")), as.list(exprs))
check(length(assign_spec) == 1L, "the script assigns the specification exactly once")
env <- new.env(); eval(assign_spec[[1]], env)
check(identical(env$spec, spec), "the embedded literal is the specification it came from")
check(isTRUE(all.equal(simulate_design(env$spec)$score, d$score)), "and reproduces the same data")

# The power curve and the sample size it solves for, the reason solve_curve.R is staged into
# the browser build at all.
cat("=== power curve and target N ===\n")
grid <- c(40, 80, 120, 160, 200)
pw <- vapply(grid, function(nn) {
  s <- spec; s$units$subject$n <- as.integer(nn); power_design(s, n_sims = 400)$power
}, numeric(1))
check(!is.unsorted(round(pw, 1)), sprintf("power rises with N (%s)", paste(round(pw, 2), collapse = ", ")))
solved <- target_n(data.frame(n_subject = grid, power = pw, n_sims = 400), target = 0.8)
check(solved$n > 90 && solved$n < 180 && solved$n_lo <= solved$n && solved$n <= solved$n_hi,
      sprintf("target_n() solves for 0.80 power (N = %d, interval %d to %d)",
              solved$n, solved$n_lo, solved$n_hi))

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

# Every response family the app offers, since each one takes a different branch of the builder
# and of the response draw. The parameters are the app's own starting values for the family.
cat("=== every family the app offers ===\n")
fam_pars <- list(
  gaussian          = list(intercept = 100, effect = 5,   sigma = 10),
  shifted_lognormal = list(intercept = 6,   effect = 0.1, sigma = 0.3, shift = 200),
  bernoulli         = list(intercept = 0,   effect = 0.5),
  poisson           = list(intercept = 1.5, effect = 0.3),
  ordinal           = list(intercept = 0,   effect = 0.8, thresholds = "-2, -0.6, 0.6, 2"),
  beta              = list(intercept = 0,   effect = 0.8, phi = 8))
for (fam in names(fam_pars)) {
  s <- build_spec(c(list(name = "ui_family", seed = 7, n_subject = 40, include_items = FALSE,
                         n_item = 24, design_kind = "between", factor_name = "group",
                         lev1 = "a", lev2 = "b", family = fam, resp_name = ""), fam_pars[[fam]]))
  dd <- simulate_design(s)
  y <- dd[[s$response$name]]
  in_range <- switch(fam,
                     bernoulli = all(y %in% c(0, 1)),
                     poisson   = all(y >= 0 & y == round(y)),
                     ordinal   = all(y %in% seq_len(length(s$response$thresholds) + 1L)),
                     beta      = all(y > 0 & y < 1),
                     TRUE)
  check(nrow(dd) == 40 && !anyNA(y) && isTRUE(in_range),
        sprintf("%s simulates 40 rows into '%s'", fam, s$response$name))
}

cat(if (ok) "\nALL APP-LOGIC CHECKS PASSED\n" else "\nSOME CHECKS FAILED\n")
quit(status = if (ok) 0 else 1)
