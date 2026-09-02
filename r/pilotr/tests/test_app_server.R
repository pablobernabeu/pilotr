# Headless test of the live Shiny reactive graph via shiny::testServer (no browser). Drives
# the real server in the installed-package app dir: sets inputs, checks the JSON output,
# triggers Simulate and the power analysis, synchronously from source and in a worker at the end.
# The serverless variant carries its own copy of the server and is driven in
# tests/test_app_lite_server.R.

library(shiny)
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
app_dir <- file.path(here, "..", "inst", "app")

ok <- TRUE
check <- function(cond, msg) { cat(if (cond) "  [PASS] " else "  [FAIL] ", msg, "\n", sep = ""); ok <<- ok && cond }

# The app carries the theme the browser variant carries. A bare fluidPage() renders under
# Shiny's default Bootstrap 3, so the same design met the reader in a different blue depending
# on which variant they reached, and the Bootstrap 5 spacing utility under the power buttons
# did nothing. The UI is built with shinyApp() stubbed out, in a fresh process: loading the app
# for its interface alone would leave the engine in this session's global environment, and the
# blocks below read that state to decide how the app finds its sources.
if (requireNamespace("bslib", quietly = TRUE) && requireNamespace("callr", quietly = TRUE)) {
  th <- tryCatch(callr::r(function(app_dir) {
    library(shiny); library(ggplot2)
    setwd(app_dir)
    env <- new.env(parent = globalenv())
    assign("shinyApp", function(ui, server) list(ui = ui, server = server), envir = env)
    built <- eval(parse("app.R"), envir = env)
    bs <- Filter(function(d) identical(d$name, "bootstrap"), htmltools::findDependencies(built$ui))
    list(bootstrap = if (length(bs) == 1L) bs[[1]]$version else NA_character_,
         primary = unname(bslib::bs_get_variables(env$APP_THEME, "primary")),
         html = as.character(htmltools::renderTags(built$ui)$html))
  }, args = list(app_dir = normalizePath(app_dir))), error = function(e) e)
  check(!inherits(th, "error") && startsWith(th$bootstrap, "5."),
        "the app renders under Bootstrap 5, as the lite app does")
  check(!inherits(th, "error") && identical(th$primary, "#2C6FB0"),
        "and on the family primary rather than Bootstrap's own")
  # A label with no text in it names its control for nobody: the paste box carried one, and a
  # screen reader was left with the placeholder, which is not a name.
  empty_labels <- if (inherits(th, "error")) "error"
    else regmatches(th$html, gregexpr("<label[^>]*></label>", th$html))[[1]]
  check(length(empty_labels) == 0L,
        "every labelled control the app renders has text in its label")
} else {
  cat("  [SKIP] theme checks (need bslib and callr)\n")
}

# A three-group between-subjects design. The point-and-click controls only ever build two
# levels, so the paste box is the only way to reach one.
three_lev <- paste0(
  '{"name":"three","seed":1,"units":{"subject":{"n":60}},',
  '"factors":[{"name":"g","levels":["a","b","c"],"contrasts":{"g1":[-1,0,1]},',
  '"between":"subject"}],"fixed":{"intercept":100,"coefficients":{"g1":5}},',
  '"random":{},"response":{"family":"gaussian","name":"score","sigma":10}}')

testServer(app = app_dir, {
  session$setInputs(
    name = "t", seed = 2024, n_subject = 64, design_kind = "between",
    include_items = FALSE, n_item = 24, factor_name = "group",
    lev1 = "control", lev2 = "treatment", intercept = 100, effect = 5,
    family = "gaussian", resp_name = "", sigma = 10)
  check(jsonlite::validate(output$json), "server renders valid JSON spec")
  # The line the Power tab opens with describes this app, which runs from the installed package
  # or from the sources beside it: it used to tell its reader to install what they already have.
  check(grepl("power_design() takes any count", output$power_out, fixed = TRUE) &&
          !grepl("install the package", output$power_out, fixed = TRUE),
        "the Power tab's opening line describes the app's own cap")
  session$setInputs(simulate = 1)
  d <- snapshot()$data
  check(nrow(d) == 64 && all(c("subject", "group", "score") %in% names(d)), "Simulate produces the 64-row data set")

  # Every tab describes the snapshot Simulate captured. Reading the response and factor names
  # from the live specification while reading the data from the last click meant that any change
  # renaming a column, a new family among them, failed with a raw R message.
  session$setInputs(family = "bernoulli")
  summ <- tryCatch(output$summary, error = function(e) paste("SUMMARY ERROR:", conditionMessage(e)))
  check(grepl("Mean (SD) of score by group", summ, fixed = TRUE),
        "the summary still describes the snapshot after the family changes")
  check(tryCatch({ output$plot; TRUE }, error = function(e) FALSE),
        "the plot still describes the snapshot after the family changes")
  check(grepl("The design has changed", output$dims, fixed = TRUE),
        "the Data tab says the design has moved on since Simulate")
  session$setInputs(family = "gaussian", simulate = 2)
  check(!grepl("The design has changed", output$dims, fixed = TRUE),
        "simulating again clears the notice")
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

  # advanced: paste a continuous-predictor spec (continuous predictors + interactions) to override
  spec_txt <- paste(readLines(file.path(here, "..", "..", "..", "spec", "examples",
                                        "reading_time_continuous.json")), collapse = "\n")
  session$setInputs(spec_json_in = spec_txt, simulate = 3)
  di <- snapshot()$data
  check("SyntaxPC" %in% names(di) && nrow(di) == 4000,
        "advanced paste-spec path simulates a continuous-predictor design")

  # An error inside an observer ends the Shiny session, taking the design with it, so the
  # power buttons report every engine refusal as text. A three-level between factor reaches
  # the guard; two subjects (one per group) reaches power_design()'s own precondition.
  session$setInputs(spec_json_in = "", n_subject = 2, run_power = 5)
  check(grepl("at least 2 subjects", output$power_out),
        "two subjects report the engine's refusal instead of ending the session")
  session$setInputs(run_curve = 1)
  check(grepl("Power curve", output$power_out), "the curve still runs from its own grid")
  check(tryCatch({ output$power_plot; TRUE }, error = function(e) FALSE), "the curve draws")

  session$setInputs(spec_json_in = three_lev, run_power = 6)
  check(grepl("two-group Gaussian design", output$power_out),
        "a three-level between factor gets the unsupported-design message")

  # The fill scales take as many colours as the data has levels. A fixed pair of them meant
  # that the same three-level design, which only the paste box can express, replaced the
  # Summary & plot tab with "Insufficient values in manual scale".
  session$setInputs(simulate = 4)
  check(tryCatch({ output$plot; TRUE }, error = function(e) FALSE),
        "a three-level pasted design plots rather than exhausting the fill scale")
  seven_lev <- paste0(
    '{"name":"seven","seed":1,"units":{"subject":{"n":70}},',
    '"factors":[{"name":"g","levels":[', paste(sprintf('"L%d"', 1:7), collapse = ","), '],',
    '"contrasts":{"g1":[-3,-2,-1,0,1,2,3]},"between":"subject"}],',
    '"fixed":{"intercept":100,"coefficients":{"g1":2}},',
    '"random":{},"response":{"family":"gaussian","name":"score","sigma":10}}')
  session$setInputs(spec_json_in = seven_lev, simulate = 5)
  check(tryCatch({ output$plot; TRUE }, error = function(e) FALSE),
        "more levels than the palette holds still plots")

  # A pasted spec that does not parse: the observers and the downloads report it, where a
  # validate() call would abort the observers silently and fail the download in the browser.
  session$setInputs(spec_json_in = "{ not json ", run_power = 7)
  check(grepl("lexical error", output$power_out), "a malformed paste is reported by Run power")
  session$setInputs(verify_code = 2)
  check(grepl("lexical error", output$verify_out), "a malformed paste is reported by Verify")
  check(identical(readLines(output$dl_spec, warn = FALSE), "{}"),
        "the spec download degrades to a placeholder rather than failing")
  session$setInputs(spec_json_in = "", n_subject = 64, run_power = 8)
  check(grepl("Power", output$power_out), "a corrected specification runs again")

  # A name that would take another column's place, and a blank one, are refused before any data
  # is made: the response used to be written over the factor column without a word, and an
  # emptied Factor name box stopped inside simulate_design() with "replacement has length zero".
  session$setInputs(resp_name = "group", simulate = 6)
  msg <- tryCatch({ output$dims; "no message" }, error = function(e) conditionMessage(e))
  check(grepl("is used by factors[1].name and response.name", msg, fixed = TRUE),
        "a response named after the factor is reported rather than simulated")
  session$setInputs(resp_name = "", factor_name = "", simulate = 7)
  msg <- tryCatch({ output$dims; "no message" }, error = function(e) conditionMessage(e))
  check(grepl("factors[1].name must be a non-empty string", msg, fixed = TRUE),
        "an emptied Factor name box is reported by the field it belongs to")
  session$setInputs(factor_name = "group", simulate = 8)

  # The numeric boxes. A cleared one arrives as NA, and as.integer() carried it into the spec as
  # NA_integer_, which jsonlite writes as the string "NA": the Design spec tab showed it and the
  # Download spec button wrote it, and only load_spec() refused it, away from the app. A typed
  # decimal was truncated in silence, so the data came from a seed nobody chose.
  session$setInputs(seed = NA_integer_)
  msg <- tryCatch({ output$json; "no message" }, error = function(e) conditionMessage(e))
  check(grepl("'seed' must be a single whole number", msg, fixed = TRUE),
        "a cleared Seed box is reported by the Design spec tab")
  check(identical(readLines(output$dl_spec, warn = FALSE), "{}"),
        "the spec download refuses with the tab rather than writing a seed of NA")
  session$setInputs(seed = 2024.5)
  msg <- tryCatch({ output$json; "no message" }, error = function(e) conditionMessage(e))
  check(grepl("'seed' must be a single whole number", msg, fixed = TRUE),
        "a fractional seed is reported rather than truncated to another design")

  # Download names come from the specification, not from the sidebar: with a pasted spec in force
  # those are two different designs, and an emptied Design name wrote '.json', '.R' and
  # '_seed2024.csv', which are hidden on Unix and refused by some browsers.
  session$setInputs(seed = 2024, name = "")
  check(identical(basename(output$dl_spec), "design.json") &&
          identical(basename(output$dl_rscript), "design.R") &&
          identical(basename(output$dl_data), "design_seed2024.csv"),
        "an emptied Design name falls back rather than naming dot-files")
  session$setInputs(name = "my_design", spec_json_in = three_lev)
  check(identical(basename(output$dl_spec), "three.json") &&
          identical(basename(output$dl_rscript), "three.R") &&
          identical(basename(output$dl_data), "three_seed1.csv"),
        "a pasted design's downloads carry its own name and seed")
  session$setInputs(spec_json_in = "", name = "t")

  # verified R-script export: run the script the tab shows in a clean R subprocess and compare
  session$setInputs(spec_json_in = "", verify_code = 1)
  vo <- output$verify_out
  check(grepl("bit for bit", vo, fixed = TRUE), "verify: the exported script reproduces the data bit for bit")
  cat("  verify:", gsub("\n", " ", vo), "\n")

  # What runs in the clean session is the script, not the specification behind it. While Verify
  # ran the JSON, a defect confined to the script emitter, which quotes list names and formats
  # its own numbers, left the button green. The app sources the engine into the global
  # environment, which is where the observer resolves the emitter from.
  gen_ok <- generate_r_script
  assign("generate_r_script", function(spec)
    sub("intercept = 100", "intercept = 101", gen_ok(spec), fixed = TRUE), envir = globalenv())
  session$setInputs(verify_code = 3)
  check(grepl("not identical", output$verify_out, fixed = TRUE),
        "verify: a defect in the script emitter is caught rather than passed over")
  assign("generate_r_script", gen_ok, envir = globalenv())
})

# Both power buttons take the same route to the background worker that run_app(async = TRUE)
# sets up. The curve is the point estimate repeated once per sample size, and it ran in the main
# process while the cheaper button went to the worker, so the session froze for the whole sweep.
# The branch is only reachable with pilotr as a loaded namespace (a worker reloads the package),
# which the blocks above deliberately avoid, so it is driven last and only where the package can
# be loaded; a sequential plan keeps the worker in this process.
if ((nzchar(system.file(package = "pilotr")) || requireNamespace("pkgload", quietly = TRUE)) &&
      requireNamespace("future", quietly = TRUE) && requireNamespace("promises", quietly = TRUE)) {
  # Shiny sources an app directory into the global environment, so the blocks above left the
  # engine there and loading the package reports it as a conflict. The app runs the same code
  # either way; what the namespace adds is the async branch.
  if (!nzchar(system.file(package = "pilotr")))
    suppressWarnings(suppressMessages(pkgload::load_all(file.path(here, ".."), quiet = TRUE,
                                                        export_all = FALSE)))
  else suppressMessages(library(pilotr))
  future::plan(future::sequential)
  testServer(app = app_dir, {
    session$setInputs(
      name = "t", seed = 2024, n_subject = 20, design_kind = "between",
      include_items = FALSE, n_item = 24, factor_name = "group",
      lev1 = "control", lev2 = "treatment", intercept = 100, effect = 5,
      family = "gaussian", resp_name = "", sigma = 10, n_sims = 100)
    session$setInputs(run_curve = 1)
    check(grepl("background worker", output$power_out),
          "async: the curve hands its sweep to the worker rather than blocking the session")
    session$flushReact()
    check(grepl("Power curve at n_sims = 100", output$power_out),
          "async: the worker's result reaches the Power tab")
    check(tryCatch({ output$power_plot; TRUE }, error = function(e) FALSE), "async: the curve draws")
    session$setInputs(run_power = 1)
    check(grepl("background worker", output$power_out), "async: so does the point estimate")
    session$flushReact()
    check(grepl("Power +:", output$power_out), "async: and its result reaches the Power tab")
  })
  future::plan(future::sequential)
} else {
  cat("  [SKIP] async power runs (needs pilotr as a namespace, future and promises)\n")
}

cat(if (ok) "TESTSERVER OK\n" else "TESTSERVER FAILED\n")
quit(status = if (ok) 0 else 1)
