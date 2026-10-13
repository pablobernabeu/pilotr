# Headless test of the browser build's server via shiny::testServer, with no browser. The build is
# staged as app-lite/build_shinylive.R stages it, with the engine files that
# app-lite/engine-files.txt names copied beside app-lite/app.R, so the app runs on the engine
# subset the site ships and on nothing else. tests/test_app_server.R makes the same checks of the
# installed app, and checks the text both apps give the designs their power tabs cannot analyse.

library(shiny)
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
lite_src <- file.path(here, "..", "..", "..", "app-lite")

ok <- TRUE
check <- function(cond, msg) { cat(if (cond) "  [PASS] " else "  [FAIL] ", msg, "\n", sep = ""); ok <<- ok && cond }

stage <- file.path(tempdir(), "pilotr-lite")
unlink(stage, recursive = TRUE)
dir.create(stage)
engine <- trimws(sub("#.*$", "", readLines(file.path(lite_src, "engine-files.txt"))))
engine <- engine[nzchar(engine)]
stopifnot(all(file.copy(file.path(here, "..", "R", engine), stage)),
          file.copy(file.path(lite_src, "app.R"), stage))

# Each check below starts a session of its own, from the values the controls start at with fewer
# simulations. An error in an observer closes the session, and a check that followed it would
# read a closed session.
start_inputs <- list(
  name = "t", seed = 2024, n_subject = 64, design_kind = "between", include_items = FALSE,
  n_item = 24, factor_name = "group", lev1 = "control", lev2 = "treatment", intercept = 100,
  effect = 5, family = "gaussian", resp_name = "", sigma = 10, n_sims = 100,
  use_pasted = FALSE, pasted = "")

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
# seed of its own, 99. build_spec() is the staged engine's, which the app sources into the global
# environment.
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
# grid point.
survival <- list(
  list(what = "an N of 2", set = list(n_subject = 2),
       says = "Set N subjects to at least 4 before running power."),
  list(what = "a cleared N", set = list(n_subject = NA),
       says = "Set N subjects to at least 4 before running power."),
  list(what = "a pasted three-level design", paste = "three",
       says = "The in-browser demo runs power only for the two-group Gaussian design"),
  list(what = "a pasted 2 x 2 design", paste = "two_by_two",
       says = "are both between 'subject'"))
buttons <- c(run_power = "Estimate power", run_curve = "Power curve")
for (case in survival) for (button in names(buttons)) {
  testServer(app = stage, {
    do.call(session$setInputs, start_inputs)
    if (!is.null(case$set)) do.call(session$setInputs, case$set)
    if (!is.null(case$paste))
      session$setInputs(use_pasted = TRUE, pasted = pasted_designs()[[case$paste]])
    do.call(session$setInputs, setNames(list(1L), button))
    check(!session$isClosed() && grepl(case$says, shown(output$power_out), fixed = TRUE),
          sprintf("browser app: %s keeps the session open when %s is selected, and says why",
                  case$what, buttons[[button]]))
  })
}

# Data, Summary & plot and the data download show the last simulation, labelled with its own seed,
# and a note says when the controls have changed since. Before, a new family or response name
# stopped the summary and the plot, and a new factor name dropped their grouping. The Data tab gave
# no seed, and the download simulated afresh from the controls.
moved_on <- paste0(
  "The controls have changed since the last Simulate: the JSON and the specification download ",
  "follow the controls, while the data, summary, plot and data download show the last ",
  "simulation. Select Simulate to bring them together.")
testServer(app = stage, {
  do.call(session$setInputs, start_inputs)
  session$setInputs(simulate = 1)
  simulated <- simulate_design(current_spec())
  check(identical(shown(output$sim_note), ""),
        "browser app: no note while the controls match the last simulation")
  changes <- list(list(family = "bernoulli"), list(resp_name = "rt"), list(factor_name = "arm"))
  for (change in changes) {
    do.call(session$setInputs, change)
    label <- sprintf("%s = \"%s\"", names(change), change[[1]])
    summary_text <- shown(output$summary)
    check(grepl("treatment", summary_text, fixed = TRUE) &&
            grepl("mean", summary_text, fixed = TRUE),
          sprintf("browser app: after %s, the summary describes the last simulation by group",
                  label))
    check(identical(x_levels(output$plot), c("control", "treatment")),
          sprintf("browser app: after %s, the plot draws the last simulation by group", label))
    do.call(session$setInputs, start_inputs[names(change)])
  }
  session$setInputs(seed = 7)
  check(grepl("(seed 2024)", shown(output$dims), fixed = TRUE),
        "browser app: the Data tab gives the seed of the last simulation")
  check(grepl(moved_on, shown(output$sim_note), fixed = TRUE),
        "browser app: a note says that the controls have changed since the last Simulate")
  download <- tryCatch(output$dl_data, error = function(e) "")
  check(identical(basename(download), "t_seed2024.csv") &&
          identical(readLines(download), csv_lines(simulated)),
        "browser app: the data download holds the last simulation and is named with its seed")
  session$setInputs(simulate = 2)
  check(grepl("(seed 7)", shown(output$dims), fixed = TRUE) &&
          identical(shown(output$sim_note), ""),
        "browser app: simulating again brings the data and the controls together")
  session$setInputs(use_pasted = TRUE, pasted = pasted_designs()$three, simulate = 3)
  check(grepl("(seed 99)", shown(output$dims), fixed = TRUE),
        "browser app: a pasted specification is shown with its own seed")
  check(identical(x_levels(output$plot), c("a", "b", "c")),
        "browser app: a pasted three-level design plots its three groups")
})

cat(if (ok) "LITE TESTSERVER OK\n" else "LITE TESTSERVER FAILED\n")
quit(status = if (ok) 0 else 1)
