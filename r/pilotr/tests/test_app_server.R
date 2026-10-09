# Headless test of the live Shiny reactive graph via shiny::testServer (no browser). Drives
# the real server in the installed-package app dir: sets inputs, checks the JSON output,
# triggers Simulate and the (synchronous, from-source) power analysis. A staged copy of the
# browser build is then driven for the text it gives the designs its power tab cannot analyse.

library(shiny)
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
app_dir <- file.path(here, "..", "inst", "app")

ok <- TRUE
check <- function(cond, msg) { cat(if (cond) "  [PASS] " else "  [FAIL] ", msg, "\n", sep = ""); ok <<- ok && cond }

# The designs whose not-supported text both apps route, built with the app's own build_spec().
# They are crossed reaction times and crossed accuracy with one within factor, a Likert and a
# proportion design between subjects, and a Gaussian design crossed with items that has no within
# factor.
route_designs <- function(build_spec) {
  rt <- build_spec(list(name = "rt", seed = 1, design_kind = "within", include_items = TRUE,
                        n_subject = 12, n_item = 8, factor_name = "cond", lev1 = "a",
                        lev2 = "b", intercept = 6, effect = 0.1, subj_int_sd = 0.1,
                        subj_slope_sd = 0, item_int_sd = 0.1, item_slope_sd = 0,
                        family = "shifted_lognormal", resp_name = "", sigma = 0.3,
                        shift = 200))
  acc <- rt
  acc$response <- list(family = "bernoulli", name = "accuracy")
  acc$fixed$intercept <- 1
  between <- function(family, ...)
    build_spec(list(name = family, seed = 1, design_kind = "between", n_subject = 40,
                    factor_name = "group", lev1 = "a", lev2 = "b", intercept = 0,
                    effect = 0.8, family = family, resp_name = "", ...))
  crossed <- between("gaussian", sigma = 1)
  crossed$units$item <- list(n = 20L)
  crossed$random <- list(subject = list(intercept_sd = 1), item = list(intercept_sd = 0.3))
  list(rt = rt, acc = acc, likert = between("ordinal", thresholds = "-1, 0, 1"),
       prop = between("beta", phi = 8), crossed = crossed)
}

# The not-supported text sends each design to the analysis that fits it. power_mixed() in R fits
# every family but ordinal and Beta, which go to the Bayesian design analysis. Python's
# power_mixed() takes one within factor crossed with items and fits a linear model on the
# response's own scale. Its lines are offered only for such a design, in a family whose
# coefficients are on that scale (gaussian, ex-Gaussian) or on the log scale it analyses. `out`
# holds each design's text, `opening` the sentence each text starts from.
check_routes <- function(out, opening, app) {
  check(all(vapply(out, grepl, logical(1), pattern = opening, fixed = TRUE)),
        sprintf("%s: every routed design opens with the app's own scope", app))
  check(grepl("power_mixed(spec", out$rt, fixed = TRUE) && grepl("Python", out$rt, fixed = TRUE),
        sprintf("%s: a crossed reaction-time design with one within factor is sent to both packages", app))
  check(grepl("power_mixed(spec", out$acc, fixed = TRUE) && !grepl("Python", out$acc, fixed = TRUE),
        sprintf("%s: a crossed accuracy design is sent to R's power_mixed() alone", app))
  for (nm in c("likert", "prop"))
    check(grepl("generate_design_analysis(", out[[nm]], fixed = TRUE) &&
            grepl("brms_bridge(", out[[nm]], fixed = TRUE) &&
            !grepl("power_mixed(", out[[nm]], fixed = TRUE),
          sprintf("%s: a '%s' design is sent to the Bayesian design analysis", app, nm))
  check(grepl("power_mixed(spec", out$crossed, fixed = TRUE) &&
          !grepl("Python", out$crossed, fixed = TRUE),
        sprintf("%s: a design with no within factor is sent to R's power_mixed() alone", app))
}

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

  routes <- route_designs(build_spec)
  out <- list()
  for (nm in names(routes)) {
    clicks <- clicks + 1
    session$setInputs(spec_json_in = spec_json(routes[[nm]]), run_power = clicks)
    out[[nm]] <- output$power_out
  }
  check_routes(out, "The in-app power backend covers", "installed app")

  # verified R-script export: run the design in a clean R subprocess and compare
  session$setInputs(spec_json_in = "", verify_code = 1)
  vo <- output$verify_out
  check(grepl("reproduces identically", vo, ignore.case = TRUE), "verify: clean R session reproduces the data bit-for-bit")
  cat("  verify:", gsub("\n", " ", vo), "\n")
})

# The browser build writes the same routing out in its own copy of the text. It is staged as
# build_shinylive.R stages it, with the engine files app-lite/engine-files.txt names beside
# app-lite/app.R, and each design is pasted in as a specification.
lite <- file.path(tempdir(), "pilotr-lite")
unlink(lite, recursive = TRUE)
dir.create(lite)
lite_src <- file.path(here, "..", "..", "..", "app-lite")
engine <- trimws(sub("#.*$", "", readLines(file.path(lite_src, "engine-files.txt"))))
engine <- engine[nzchar(engine)]
stopifnot(all(file.copy(file.path(here, "..", "R", engine), lite)),
          file.copy(file.path(lite_src, "app.R"), lite))
testServer(app = lite, {
  routes <- route_designs(build_spec)
  session$setInputs(n_sims = 100, use_pasted = TRUE)
  out <- list()
  for (i in seq_along(routes)) {
    session$setInputs(pasted = spec_json(routes[[i]]), run_power = i)
    out[[names(routes)[i]]] <- output$power_out
  }
  check_routes(out, "runs power only for the two-group Gaussian design", "browser app")
})

cat(if (ok) "TESTSERVER OK\n" else "TESTSERVER FAILED\n")
quit(status = if (ok) 0 else 1)
