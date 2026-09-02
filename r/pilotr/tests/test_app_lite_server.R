# Headless test of the serverless (shinylive / webR) app's reactive graph via shiny::testServer.
#
# app-lite carries its own copy of the server, and it is the variant most people meet: the site
# workflow publishes it to GitHub Pages and the vignettes link to it as the way that needs nothing
# installed. Driving only the bundled app left every one of its 500 lines untested, so a fix
# landed in one copy and not the other; the three-level fill scale was found that way.
#
# The app is staged exactly as app-lite/build_shinylive.R stages it, from the same manifest, so
# what runs here is the build the site ships rather than the source tree.

library(shiny)
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
lite_src <- file.path(here, "..", "..", "..", "app-lite")

ok <- TRUE
check <- function(cond, msg) { cat(if (cond) "  [PASS] " else "  [FAIL] ", msg, "\n", sep = ""); ok <<- ok && cond }

stage <- file.path(tempdir(), "app-lite-staged")
unlink(stage, recursive = TRUE); dir.create(stage, showWarnings = FALSE)
engine <- trimws(sub("#.*$", "", readLines(file.path(lite_src, "engine-files.txt"))))
engine <- engine[nzchar(engine)]
paths <- file.path(here, "..", "R", engine)
stopifnot(length(paths) > 0, all(file.exists(paths)))
stopifnot(all(file.copy(paths, stage)), file.copy(file.path(lite_src, "app.R"), stage))
old_wd <- setwd(stage)   # app.R sources whatever was staged beside it
on.exit(setwd(old_wd), add = TRUE)

# A three-group between-subjects design. The point-and-click controls only ever build two levels,
# so the paste box is the only way to reach one.
three_lev <- paste0(
  '{"name":"three","seed":1,"units":{"subject":{"n":60}},',
  '"factors":[{"name":"g","levels":["a","b","c"],"contrasts":{"g1":[-1,0,1]},',
  '"between":"subject"}],"fixed":{"intercept":100,"coefficients":{"g1":5}},',
  '"random":{},"response":{"family":"gaussian","name":"score","sigma":10}}')

testServer(app = stage, {
  session$setInputs(
    name = "t", seed = 2024, n_subject = 64, design_kind = "between",
    include_items = FALSE, n_item = 24, factor_name = "group",
    lev1 = "control", lev2 = "treatment", intercept = 100, effect = 5,
    family = "gaussian", resp_name = "", sigma = 10,
    use_pasted = FALSE, pasted = "", n_sims = 300)
  check(jsonlite::validate(output$json), "the server renders valid JSON")
  session$setInputs(simulate = 1)
  d <- snapshot()$data
  check(nrow(d) == 64 && all(c("subject", "group", "score") %in% names(d)),
        "Simulate produces the 64-row data set")
  check(grepl("Simulated 64 rows x 3 columns", output$dims, fixed = TRUE), "the Data tab counts them")
  check(tryCatch({ output$plot; output$summary; TRUE }, error = function(e) FALSE),
        "the summary and the plot describe the snapshot")

  session$setInputs(run_power = 1)
  check(grepl("Power: 0", output$power_out) && grepl("n_sims = 300", output$power_out),
        "the point estimate runs at the requested count")
  session$setInputs(run_curve = 1)
  # The in-browser sweep caps its replicates, since it has no worker to run in.
  check(grepl("Power curve at n_sims = 200", output$power_out, fixed = TRUE) &&
          grepl("Target power 0.80", output$power_out, fixed = TRUE),
        "the curve runs capped and solves for the target")
  check(tryCatch({ output$power_plot; TRUE }, error = function(e) FALSE), "the curve draws")

  # A cleared Simulations box arrives as NA, and clamping NA leaves NA, which reaches
  # power_design() as 'vector size cannot be NA' and says nothing about the box that was emptied.
  session$setInputs(n_sims = NA_integer_, run_power = 2)
  check(grepl("n_sims = 300", output$power_out, fixed = TRUE),
        "a cleared Simulations box falls back to the count the box starts at")
  session$setInputs(n_sims = 5L, run_power = 3)
  check(grepl("n_sims = 100", output$power_out, fixed = TRUE),
        "a count below the minimum is clamped")
  session$setInputs(n_sims = 300)

  # The downloads are named after the specification they carry, and an emptied Design name box
  # falls back rather than writing a dot-file.
  check(identical(basename(output$dl_spec), "t.json") &&
          identical(basename(output$dl_data), "t.csv") &&
          identical(basename(output$dl_script), "t.R"),
        "the downloads carry the design's name")
  session$setInputs(name = "")
  check(identical(basename(output$dl_spec), "design.json"),
        "an emptied Design name falls back rather than naming a dot-file")
  session$setInputs(name = "t")

  # The fill scales take as many colours as the data has levels, and the power guard asks for
  # exactly one two-level between factor, so a three-level design reaches neither ggplot2's
  # 'Insufficient values in manual scale' nor the engine's own refusal.
  session$setInputs(use_pasted = TRUE, pasted = three_lev, simulate = 2)
  check(nrow(snapshot()$data) == 60, "a three-level pasted design simulates")
  check(tryCatch({ output$plot; TRUE }, error = function(e) FALSE),
        "a three-level pasted design plots rather than exhausting the fill scale")
  session$setInputs(run_power = 4)
  check(grepl("two-group Gaussian design", output$power_out),
        "a three-level between factor gets the unsupported-design message")

  # A pasted specification that does not parse: the observers and the downloads report it, where
  # an error inside an observer would end the session and a failing download would leave the
  # browser with nothing at all.
  session$setInputs(pasted = "{ not json ", run_power = 5)
  check(grepl("lexical error", output$power_out), "a malformed paste is reported by the power button")
  check(identical(readLines(output$dl_spec, warn = FALSE), "{}"),
        "the spec download degrades to a placeholder rather than failing")
  check(identical(readLines(output$dl_data, warn = FALSE), ""),
        "so does the data download")
  msg <- tryCatch({ output$json; "no message" }, error = function(e) conditionMessage(e))
  check(grepl("lexical error", msg), "and the Design spec tab says why")

  # The controls can also build a specification the package refuses: a cleared Seed box arrives
  # as NA, which writes as the string "NA" in a file load_spec() then will not read.
  session$setInputs(use_pasted = FALSE, seed = NA_integer_)
  msg <- tryCatch({ output$json; "no message" }, error = function(e) conditionMessage(e))
  check(grepl("'seed' must be a single whole number", msg, fixed = TRUE),
        "a cleared Seed box is reported by the Design spec tab")
  check(identical(readLines(output$dl_spec, warn = FALSE), "{}"),
        "the spec download refuses rather than writing a seed of NA")
  session$setInputs(seed = 2024, run_power = 6)
  check(grepl("Power: 0", output$power_out), "a corrected specification runs again")
})
setwd(old_wd)

cat(if (ok) "LITE TESTSERVER OK\n" else "LITE TESTSERVER FAILED\n")
quit(status = if (ok) 0 else 1)
