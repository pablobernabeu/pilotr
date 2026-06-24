# Headless test of the live Shiny reactive graph via shiny::testServer (no browser). Drives
# the real server in the installed-package app dir: sets inputs, checks the JSON output,
# triggers Simulate and the (synchronous, from-source) power analysis.

library(shiny)
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
app_dir <- file.path(here, "..", "inst", "app")

ok <- TRUE
check <- function(cond, msg) { cat(if (cond) "  [PASS] " else "  [FAIL] ", msg, "\n", sep = ""); ok <<- ok && cond }

testServer(app = app_dir, {
  session$setInputs(
    name = "t", seed = 2024, n_subject = 64, design_kind = "between",
    include_items = FALSE, n_item = 24, factor_name = "group",
    lev1 = "control", lev2 = "treatment", intercept = 100, effect = 5,
    family = "gaussian", resp_name = "", sigma = 10)
  check(jsonlite::validate(output$json), "server renders valid JSON spec")
  session$setInputs(simulate = 1)
  d <- data()
  check(nrow(d) == 64 && all(c("subject", "group", "score") %in% names(d)), "Simulate produces the 64-row data set")
  session$setInputs(n_sims = 300, run_power = 1)
  po <- output$power_out
  check(grepl("Power", po), "power analysis output rendered")
  cat("  power output:\n", gsub("\n", "\n    ", po), "\n")

  # advanced: paste a continuous-predictor spec (continuous predictors + interactions) to override
  spec_txt <- paste(readLines(file.path(here, "..", "..", "..", "spec", "examples",
                                        "reading_time_continuous.json")), collapse = "\n")
  session$setInputs(spec_json_in = spec_txt, simulate = 2)
  di <- data()
  check("SyntaxPC" %in% names(di) && nrow(di) == 4000,
        "advanced paste-spec path simulates a continuous-predictor design")

  # verified R-script export: run the design in a clean R subprocess and compare
  session$setInputs(spec_json_in = "", verify_code = 1)
  vo <- output$verify_out
  check(grepl("reproduces identically", vo, ignore.case = TRUE), "verify: clean R session reproduces the data bit-for-bit")
  cat("  verify:", gsub("\n", " ", vo), "\n")
})

cat(if (ok) "TESTSERVER OK\n" else "TESTSERVER FAILED\n")
quit(status = if (ok) 0 else 1)
