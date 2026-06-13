# Headless test of the live Shiny reactive graph via shiny::testServer (no browser).
# Drives the real server: sets inputs, checks the JSON output, triggers Simulate and the
# power analysis, and asserts the reactive outputs are correct.

library(shiny)
ok <- TRUE
check <- function(cond, msg) { cat(if (cond) "  [PASS] " else "  [FAIL] ", msg, "\n", sep = ""); ok <<- ok && cond }

testServer(app = file.path("toolkit", "app"), {
  session$setInputs(
    name = "t", seed = 2024, n_subject = 64, design_kind = "between",
    include_items = FALSE, n_item = 24, factor_name = "group",
    lev1 = "control", lev2 = "treatment", intercept = 100, effect = 5,
    family = "gaussian", resp_name = "", sigma = 10)

  check(jsonlite::validate(output$json), "server renders valid JSON spec")

  session$setInputs(simulate = 1)
  d <- data()
  check(nrow(d) == 64 && all(c("subject", "group", "score") %in% names(d)),
        "Simulate produces the 64-row data set")

  session$setInputs(n_sims = 300, run_power = 1)
  po <- output$power_out
  check(grepl("Power", po), "power analysis output rendered")
  cat("  power output:\n", gsub("\n", "\n    ", po), "\n")
}, args = list())

cat(if (ok) "TESTSERVER OK\n" else "TESTSERVER FAILED\n")
quit(status = if (ok) 0 else 1)
