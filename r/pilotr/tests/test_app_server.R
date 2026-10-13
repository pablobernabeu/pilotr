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
# The sign both apps print between Type M and its standard error, built here so that this file
# stays ASCII.
plus_minus <- intToUtf8(0xB1)

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

# The power at each point of an app's curve, worked out as sweep_spec() would: point j of the grid
# takes replicate_seeds(seed, J)[j] as its seed. Both apps put the engine in the global
# environment, which is where this finds power_design() and replicate_seeds().
curve_at_grid_seeds <- function(spec, grid, n_sims) {
  seeds <- replicate_seeds(spec$seed, length(grid))
  vapply(seq_along(grid), function(j) {
    s <- spec
    s$units$subject$n <- as.integer(grid[j])
    s$seed <- seeds[j]
    power_design(s, n_sims = n_sims)$power
  }, numeric(1))
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
  d <- sim()$data
  check(nrow(d) == 64 && all(c("subject", "group", "score") %in% names(d)), "Simulate produces the 64-row data set")
  session$setInputs(n_sims = 300, run_power = 1)
  po <- output$power_out
  check(grepl("Power", po), "power analysis output rendered")
  cat("  power output:\n", gsub("\n", "\n    ", po), "\n")
  # Every rate is shown with its Monte Carlo error: power with its standard error and Wilson
  # interval, Type S with its Wilson interval and Type M with its standard error. The output
  # gave the bare rates, and Type M rests on the significant replicates alone.
  r <- power_result()
  check(grepl(sprintf("Power        : %.3f (MCSE %.3f; 95%% Wilson [%.3f, %.3f])",
                      r$power, r$power_mcse, r$power_lo, r$power_hi), po, fixed = TRUE),
        "power is shown with its Monte Carlo standard error and Wilson interval")
  check(grepl(sprintf("Type S error : %.4f (95%% Wilson [%.4f, %.4f])",
                      r$type_s, r$type_s_lo, r$type_s_hi), po, fixed = TRUE),
        "Type S is shown with its Wilson interval")
  check(grepl(sprintf("Type M (exag): %.3f %s %.3f (MCSE)",
                      r$type_m, plus_minus, r$type_m_mcse), po, fixed = TRUE),
        "Type M is shown with its Monte Carlo standard error")
  # With no true effect, Type S and Type M are undefined, and so is their uncertainty, so both
  # are shown bare, as before. The click count is set apart from the ones below, since a click
  # that repeats the previous value fires nothing.
  session$setInputs(effect = 0, run_power = 101)
  po0 <- output$power_out
  check(grepl("Type S error : NaN\n", po0, fixed = TRUE) &&
          grepl("Type M (exag): NaN\n", po0, fixed = TRUE) &&
          grepl("(MCSE ", po0, fixed = TRUE),
        "a design with no true effect shows Type S and Type M bare and power with its error")
  session$setInputs(effect = 5)

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

  # The power curve gives each grid point its own seed, by sweep_spec()'s rule, so the points it
  # hands target_n() are independent. Every point used to keep the design's seed.
  session$setInputs(n_sims = 100, run_curve = 1)
  pc <- power_curve_data()
  check(identical(pc$pw, curve_at_grid_seeds(current_spec(), pc$grid, 100)),
        "the power curve gives each grid point its own seed")
  session$setInputs(n_sims = 300)

  # build_spec() keeps a seed of 2^31 or more, which it used to turn into NA. The caption
  # formatted the seed with %d, which stops on a double outside the integer range.
  session$setInputs(seed = 3e9, simulate = 2)
  check(identical(nrow(sim()$data), 64L) && grepl("(seed 3000000000)", output$dims, fixed = TRUE),
        "a seed of 3e9 simulates and the caption prints it in full")
  session$setInputs(seed = 2024, simulate = 3)

  # advanced: paste a continuous-predictor spec (continuous predictors + interactions) to override
  spec_txt <- paste(readLines(file.path(here, "..", "..", "..", "spec", "examples",
                                        "reading_time_continuous.json")), collapse = "\n")
  session$setInputs(spec_json_in = spec_txt, simulate = 4)
  di <- sim()$data
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

  # Verify runs the script the R script tab shows, in a clean R session, and compares its data,
  # column by column, with what the current specification produces. It used to simulate from
  # the JSON in that session and compare the sums of the response column with a tolerance, so it
  # certified a script it never ran.
  session$setInputs(spec_json_in = "", verify_code = 1)
  vo <- output$verify_out
  cat("  verify:", gsub("\n", " ", vo), "\n")
  check(identical(vo, sprintf(paste0(
          "Reproduces bit for bit in a clean R session: the script's data equal the data the ",
          "current specification produces (%d rows, all columns identical)."),
          nrow(simulate_design(current_spec())))),
        "verify: the exported script reproduces the data bit for bit in a clean R session")

  # A script that does not reproduce the data is caught, and the verdict names where the data
  # first differ. The app sourced the engine into the global environment, which is where Verify
  # finds generate_r_script(), so each stub is put there and the original restored afterwards.
  # The first stub changes one coefficient of the script's specification. The second reverses the
  # script's rows, which leaves unchanged the sum of the response column, all the old check
  # compared.
  real_script <- generate_r_script
  stub <- function(f) assign("generate_r_script", f, envir = globalenv())
  verdicts <- tryCatch({
    stub(function(spec) {
      spec$fixed$coefficients$effect <- spec$fixed$coefficients$effect + 1
      real_script(spec)
    })
    session$setInputs(verify_code = 2)
    coef <- output$verify_out
    stub(function(spec) paste0(real_script(spec),
                               "data <- data[rev(seq_len(nrow(data))), ]\n",
                               "rownames(data) <- NULL\n"))
    session$setInputs(verify_code = 3)
    list(coef = coef, rows = output$verify_out)
  }, finally = stub(real_script))
  cat("  verify, one coefficient changed:", verdicts$coef, "\n")
  cat("  verify, rows reversed:", verdicts$rows, "\n")
  mismatch <- paste0("Does not reproduce in a clean R session: the script's data differ from the ",
                     "data the current specification produces. The first difference is in ")
  check(identical(verdicts$coef, paste0(mismatch, "column 'score', at row 1.")),
        "verify: a script whose specification changes one coefficient is named at its response column")
  check(identical(verdicts$rows, paste0(mismatch, "column 'subject', at row 1.")),
        "verify: a script whose rows come back in reverse order is caught")

  # A specification the engine refuses is reported in the Verify panel. The reference was
  # simulated unguarded, so a cleared N subjects box ended the session at the click.
  session$setInputs(n_subject = NA, verify_code = 4)
  vo_na <- output$verify_out
  cat("  verify, N cleared:", gsub("\n", " ", vo_na), "\n")
  check(!session$isClosed() && startsWith(vo_na, "error: ") &&
          grepl("'units.subject.n' must be a whole number", vo_na, fixed = TRUE),
        "verify: a specification the engine refuses is reported, and the session stays open")
})

# The checks below start a session of their own each time, from the values the controls start at
# with fewer simulations. An error in an observer closes the session, and a check that followed
# it would read a closed session. tests/test_app_lite_server.R makes the same checks of the
# browser build.
start_inputs <- list(
  name = "t", seed = 2024, n_subject = 64, design_kind = "between", include_items = FALSE,
  n_item = 24, factor_name = "group", lev1 = "control", lev2 = "treatment", intercept = 100,
  effect = 5, family = "gaussian", resp_name = "", sigma = 10, n_sims = 100, spec_json_in = "")

# An output as text, or the error it stops with.
shown <- function(value)
  tryCatch(paste(unlist(value), collapse = ""),
           error = function(e) paste("ERROR:", conditionMessage(e)))
# The levels a plot draws along its x axis, from the coordinate map Shiny keeps for clicks and
# brushes, or NULL when the plot stops. A plot by group puts the factor's levels there. After a
# change of factor name, the plot used to draw the response there, without its grouping.
x_levels <- function(value)
  tryCatch(as.character(unlist(value$coordmap$panels[[1]]$domain$discrete_limits$x)),
           error = function(e) NULL)
# The lines write.csv() gives a data frame, which is how the app writes the data download.
csv_lines <- function(d) {
  tf <- tempfile(fileext = ".csv")
  on.exit(unlink(tf))
  write.csv(d, tf, row.names = FALSE)
  readLines(tf)
}

# Two Gaussian designs between subjects, as JSON to paste, that the power tab cannot analyse. A
# factor with three levels is refused by power_design(), and two factors with two levels each are
# refused by validate_spec(), since pilotr would confound them. The three-level design carries a
# seed of its own, 99. build_spec() is the engine's, which the app put in the global environment.
pasted_designs <- function() {
  base <- build_spec(list(name = "pasted", seed = 99, design_kind = "between", n_subject = 60,
                          factor_name = "group", lev1 = "a", lev2 = "b", intercept = 0,
                          effect = 0.5, family = "gaussian", resp_name = "y", sigma = 1))
  three <- base
  three$factors[[1]]$levels <- c("a", "b", "c")
  three$factors[[1]]$contrasts <- list(effect = c(-1, 0, 1))
  two_by_two <- base
  two_by_two$factors[[2]] <- list(name = "dose", levels = c("low", "high"),
                                  contrasts = list(dose = c(-0.5, 0.5)), between = "subject")
  list(three = spec_json(three), two_by_two = spec_json(two_by_two))
}

# The power buttons keep the session open and say why they cannot run. An N of 2 stopped
# power_design() in t.test() with "not enough observations", and a cleared N reached it as NA.
# Either error ended the session, and with it the design. On the curve, an N of 2 left a single
# grid point. A pasted design that load_spec() refuses did nothing, since the refusal reached the
# observer as a validation error, which Shiny drops without a word.
survival <- list(
  list(what = "an N of 2", set = list(n_subject = 2),
       says = "Set N subjects to at least 4 before running power."),
  list(what = "a cleared N", set = list(n_subject = NA),
       says = "Set N subjects to at least 4 before running power."),
  list(what = "a pasted three-level design", paste = "three",
       says = "The in-app power backend covers"),
  list(what = "a pasted 2 x 2 design", paste = "two_by_two",
       says = "are both between 'subject'"))
buttons <- c(run_power = "Run power", run_curve = "Power curve")
for (case in survival) for (button in names(buttons)) {
  testServer(app = app_dir, {
    do.call(session$setInputs, start_inputs)
    if (!is.null(case$set)) do.call(session$setInputs, case$set)
    if (!is.null(case$paste)) session$setInputs(spec_json_in = pasted_designs()[[case$paste]])
    do.call(session$setInputs, setNames(list(1L), button))
    check(!session$isClosed() && grepl(case$says, shown(output$power_out), fixed = TRUE),
          sprintf("installed app: %s keeps the session open when %s is selected, and says why",
                  case$what, buttons[[button]]))
  })
}

# Data, Summary & plot and the data download show the last simulation, labelled with its own seed,
# and a note says when the controls have changed since. Before, a new family or response name
# stopped the summary and the plot, and a new factor name dropped their grouping. The Data tab
# printed the Seed box, and the download simulated afresh from the controls.
moved_on <- paste0(
  "The controls have changed since the last Simulate: the JSON and the specification download ",
  "follow the controls, while the data, summary, plot and data download show the last ",
  "simulation. Select Simulate to bring them together.")
testServer(app = app_dir, {
  do.call(session$setInputs, start_inputs)
  session$setInputs(simulate = 1)
  simulated <- simulate_design(current_spec())
  check(identical(shown(output$sim_note), ""),
        "installed app: no note while the controls match the last simulation")
  changes <- list(list(family = "bernoulli"), list(resp_name = "rt"), list(factor_name = "arm"))
  for (change in changes) {
    do.call(session$setInputs, change)
    label <- sprintf("%s = \"%s\"", names(change), change[[1]])
    check(grepl("Mean (SD) of score by group", shown(output$summary), fixed = TRUE),
          sprintf("installed app: after %s, the summary describes the last simulation", label))
    check(identical(x_levels(output$plot), c("control", "treatment")),
          sprintf("installed app: after %s, the plot draws the last simulation by group", label))
    do.call(session$setInputs, start_inputs[names(change)])
  }
  session$setInputs(seed = 7)
  check(grepl("(seed 2024)", shown(output$dims), fixed = TRUE),
        "installed app: the Data tab gives the seed of the last simulation, not the Seed box's")
  check(grepl(moved_on, shown(output$sim_note), fixed = TRUE),
        "installed app: a note says that the controls have changed since the last Simulate")
  download <- tryCatch(output$dl_data, error = function(e) "")
  check(identical(basename(download), "t_seed2024.csv") &&
          identical(readLines(download), csv_lines(simulated)),
        "installed app: the data download holds the last simulation and is named with its seed")
  session$setInputs(simulate = 2)
  check(grepl("(seed 7)", shown(output$dims), fixed = TRUE) &&
          identical(shown(output$sim_note), ""),
        "installed app: simulating again brings the data and the controls together")
  session$setInputs(spec_json_in = pasted_designs()$three, simulate = 3)
  check(grepl("(seed 99)", shown(output$dims), fixed = TRUE),
        "installed app: a pasted specification is shown with its own seed")
  check(identical(x_levels(output$plot), c("a", "b", "c")),
        "installed app: a pasted three-level design plots its three groups")
})

# Counterbalanced lists, offered for a within design with items, write the specification 0.4 field.
# A box that a change of design hides keeps its value, so build_spec() ignores it where it cannot
# apply. tests/test_app_lite_server.R makes the same checks of the browser build.
within_inputs <- list(
  design_kind = "within", include_items = TRUE, n_subject = 12, n_item = 8, intercept = 6,
  effect = 0.1, sigma = 0.3, subj_int_sd = 0.1, subj_slope_sd = 0.05, subj_corr = 0,
  item_int_sd = 0.1, item_slope_sd = 0.05, item_corr = 0)
testServer(app = app_dir, {
  do.call(session$setInputs, start_inputs)
  do.call(session$setInputs, within_inputs)
  session$setInputs(counterbalance = TRUE, simulate = 1)
  json <- shown(output$json)
  check(grepl("\"counterbalance\": true", json, fixed = TRUE) &&
          grepl("\"spec_version\": \"0.4\"", json, fixed = TRUE),
        "installed app: Counterbalanced lists writes the field and declares 0.4")
  d <- sim()$data
  check(nrow(d) == 12 * 8 && anyDuplicated(d[c("subject", "item")]) == 0,
        "installed app: a counterbalanced design has one row per subject-item pair")
  session$setInputs(design_kind = "between")
  json <- shown(output$json)
  check(!grepl("counterbalance", json, fixed = TRUE) &&
          grepl("\"spec_version\": \"0.3\"", json, fixed = TRUE),
        "installed app: a between design leaves out the hidden, ticked Counterbalanced lists")
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

  # A design the browser build can analyse shows every rate with its Monte Carlo error, as the
  # installed app does. Both apps source their engine into the global environment, so the full
  # engine the installed app loaded above is still there. test_app_logic.R is what holds the
  # staged subset to the files power_design() needs.
  two <- build_spec(list(name = "two", seed = 2024, design_kind = "between", n_subject = 64,
                         factor_name = "group", lev1 = "control", lev2 = "treatment",
                         intercept = 100, effect = 5, family = "gaussian", resp_name = "",
                         sigma = 10))
  po <- tryCatch({
    session$setInputs(pasted = spec_json(two), run_power = length(routes) + 1L)
    output$power_out
  }, error = function(e) paste("error:", conditionMessage(e)))
  cat("  browser power output:\n", gsub("\n", "\n    ", po), "\n")
  r <- power_design(two, n_sims = 100)
  check(grepl(sprintf("Power: %.3f (MCSE %.3f; 95%% Wilson [%.3f, %.3f])",
                      r$power, r$power_mcse, r$power_lo, r$power_hi), po, fixed = TRUE),
        "browser app: power is shown with its Monte Carlo standard error and Wilson interval")
  check(grepl(sprintf("Type S: %.4f (95%% Wilson [%.4f, %.4f])",
                      r$type_s, r$type_s_lo, r$type_s_hi), po, fixed = TRUE),
        "browser app: Type S is shown with its Wilson interval")
  check(grepl(sprintf("Type M: %.3f %s %.3f (MCSE)", r$type_m, plus_minus, r$type_m_mcse),
              po, fixed = TRUE),
        "browser app: Type M is shown with its Monte Carlo standard error")
  session$setInputs(run_curve = 1)
  pc <- power_plot()
  check(identical(pc$pw, curve_at_grid_seeds(current_spec(), pc$grid, 100)),
        "browser app: the power curve gives each grid point its own seed")
})

cat(if (ok) "TESTSERVER OK\n" else "TESTSERVER FAILED\n")
quit(status = if (ok) 0 else 1)
