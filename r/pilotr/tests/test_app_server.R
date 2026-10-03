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

  # A cleared Simulations box arrives as NA. Clamping it left NA, and power_design() then
  # stopped with "vector size cannot be NA", which says nothing about the emptied box. The
  # count now falls back to the value the box starts at, 1000.
  session$setInputs(n_sims = NA_integer_, run_power = 2)
  po_na <- output$power_out
  check(grepl("Power", po_na), "a cleared Simulations box still runs a power analysis")
  check(grepl("Simulations *: *1000", po_na),
        "a cleared Simulations box falls back to the box's own default of 1000")

  # The same for a value outside the app's range: clamped, never passed through.
  session$setInputs(n_sims = 5L, run_power = 3)
  check(grepl("Simulations *: *100\\b", output$power_out),
        "a count below the minimum is clamped to 100")
  session$setInputs(n_sims = 300, run_power = 4)

  # build_spec() keeps a seed of 2^31 or more, which it used to turn into NA. The caption
  # formatted the seed with %d, which stops on a double outside the integer range.
  session$setInputs(seed = 3e9, simulate = 2)
  check(identical(nrow(data()), 64L) && grepl("(seed 3000000000)", output$dims, fixed = TRUE),
        "a seed of 3e9 simulates and the caption prints it in full")
  session$setInputs(seed = 2024, simulate = 3)

  # advanced: paste a continuous-predictor spec (continuous predictors + interactions) to override
  spec_txt <- paste(readLines(file.path(here, "..", "..", "..", "spec", "examples",
                                        "reading_time_continuous.json")), collapse = "\n")
  session$setInputs(spec_json_in = spec_txt, simulate = 4)
  di <- data()
  check("SyntaxPC" %in% names(di) && nrow(di) == 4000,
        "advanced paste-spec path simulates a continuous-predictor design")

  # Power runs only where the t-test holds, with one row per subject and no shared cluster. The
  # app writes power_design()'s rule out, since the installed app reaches only pilotr's exports,
  # so the two are held together here. Testing the first factor alone let a pasted crossed design
  # through to a t-test of its correlated rows, and let a three-level design stop the observer.
  base <- build_spec(list(name = "b", seed = 1, design_kind = "between", n_subject = 30,
                          factor_name = "group", lev1 = "a", lev2 = "b", intercept = 0,
                          effect = 0.5, family = "gaussian", resp_name = "y", sigma = 1))
  crossed <- base
  crossed$units$item <- list(n = 20L)
  crossed$random <- list(subject = list(intercept_sd = 1), item = list(intercept_sd = 0.3))
  sites <- base
  sites$random <- list(site = list(intercept_sd = 0.5, over = "subject", n = 6L))
  within <- base
  within$factors[[2]] <- list(name = "block", levels = c("x", "y"),
                              contrasts = list(blk = c(-0.5, 0.5)), vary_within = "subject")
  three <- base
  three$factors[[1]]$levels <- c("a", "b", "c")
  three$factors[[1]]$contrasts <- list(effect = c(-1, 0, 1))
  two_by_two <- base
  two_by_two$factors[[2]] <- list(name = "dose", levels = c("low", "high"),
                                  contrasts = list(dose = c(-0.5, 0.5)), between = "subject")
  by_subject <- base
  by_subject$random <- list(subject = list(intercept_sd = 0.5))
  # The last three reach the clauses the others miss. `lognormal` changes the family, and `both`
  # makes the factor vary within subjects as well. `item_entry` adds an item entry, which groups
  # nothing in a design without an item unit.
  lognormal <- base
  lognormal$response <- list(family = "lognormal", name = "y", sigma = 0.3)
  both <- base
  both$factors[[1]]$vary_within <- "subject"
  item_entry <- base
  item_entry$random <- list(item = list(intercept_sd = 0.3))
  designs <- list(base = base, by_subject = by_subject, crossed = crossed, sites = sites,
                  within = within, three = three, two_by_two = two_by_two,
                  lognormal = lognormal, both = both, item_entry = item_entry)
  agree <- vapply(designs, function(s)
    identical(gaussian_two_group(s), is.null(.two_group_refusal(s))), logical(1))
  check(all(agree), paste("the app's power check agrees with power_design()'s on",
                          paste(names(designs), collapse = ", ")))
  # Each refused design follows a power result, so a message left over from an observer that
  # stopped cannot pass for the refusal.
  session$setInputs(n_sims = 100)
  clicks <- 4
  for (nm in c("three", "by_subject", "crossed")) {
    clicks <- clicks + 1
    session$setInputs(spec_json_in = spec_json(designs[[nm]]), run_power = clicks)
    if (identical(nm, "by_subject"))
      check(grepl("Power *: ", output$power_out),
            "a pasted design with one row per subject and a by-subject intercept runs")
    else
      check(grepl("The in-app power backend covers", output$power_out, fixed = TRUE),
            sprintf("a pasted '%s' design gets the not-supported text", nm))
  }

  # verified R-script export: run the design in a clean R subprocess and compare
  session$setInputs(spec_json_in = "", verify_code = 1)
  vo <- output$verify_out
  check(grepl("reproduces identically", vo, ignore.case = TRUE), "verify: clean R session reproduces the data bit-for-bit")
  cat("  verify:", gsub("\n", " ", vo), "\n")
})

cat(if (ok) "TESTSERVER OK\n" else "TESTSERVER FAILED\n")
quit(status = if (ok) 0 else 1)
