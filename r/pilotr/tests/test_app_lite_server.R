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
# The power curve starts from half, three quarters, one, one and a half and twice N subjects. Each
# point takes the next of replicate_seeds(seed, 10) in the order the points are computed, which is
# the starting grid in increasing order, then any points added above it in increasing order, then
# any added below it in decreasing order. These are the seeds the points of `grid` should carry by
# that rule. replicate_seeds() is the staged engine's, which the app sources into the global
# environment.
expected_seeds <- function(seed, n, grid) {
  start <- unique(pmax(4, round(n * c(0.5, 0.75, 1, 1.5, 2))))
  run_order <- c(start, sort(grid[grid > max(start)]),
                 sort(grid[grid < min(start)], decreasing = TRUE))
  replicate_seeds(seed, 10)[match(grid, run_order)]
}

# Stand-ins for power_design() of known power at n subjects, as in tests/test_app_server.R, which
# says why none is flat. `crossing` is the normal approximation to the power of the default design.
stand_ins <- list(
  below = function(n) 0.1 + n / 10000,
  above = function(n) 0.95 + n / 10000,
  crossing = function(n) pnorm(0.25 * sqrt(n) - 1.96))

# Puts a stand-in for power_design() in the global environment, where the app sourced the engine,
# and records the size and seed of each call in `calls`. The app sources the engine there again
# whenever a session starts, so this is called inside the session, and `restore` returns the
# engine's own function.
stub_power <- function(power, calls) {
  real <- get("power_design", envir = globalenv())
  calls$n <- numeric(0)
  calls$seed <- numeric(0)
  assign("power_design", function(spec, n_sims = 1000, ...) {
    calls$n <- c(calls$n, spec$units$subject$n)
    calls$seed <- c(calls$seed, spec$seed)
    list(power = power(spec$units$subject$n))
  }, envir = globalenv())
  function() assign("power_design", real, envir = globalenv())
}

# The analytic sample size of the default design, which the app prints beside the curve's solve.
analytic_default <- 2 * stats::power.t.test(delta = 5, sd = 10, power = 0.8)$n

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

# Counterbalanced lists, offered for a within design with items, write the specification 0.4 field.
# A box that a change of design hides keeps its value, so build_spec() ignores it where it cannot
# apply. tests/test_app_server.R makes the same checks of the installed app.
within_inputs <- list(
  design_kind = "within", include_items = TRUE, n_subject = 12, n_item = 8, intercept = 6,
  effect = 0.1, sigma = 0.3, subj_int_sd = 0.1, subj_slope_sd = 0.05, subj_corr = 0,
  item_int_sd = 0.1, item_slope_sd = 0.05, item_corr = 0)
testServer(app = stage, {
  do.call(session$setInputs, start_inputs)
  do.call(session$setInputs, within_inputs)
  session$setInputs(counterbalance = TRUE, simulate = 1)
  json <- shown(output$json)
  check(grepl("\"counterbalance\": true", json, fixed = TRUE) &&
          grepl("\"spec_version\": \"0.4\"", json, fixed = TRUE),
        "browser app: Counterbalanced lists writes the field and declares 0.4")
  d <- sim()$data
  check(nrow(d) == 12 * 8 && anyDuplicated(d[c("subject", "item")]) == 0,
        "browser app: a counterbalanced design has one row per subject-item pair")
  session$setInputs(design_kind = "between")
  json <- shown(output$json)
  check(!grepl("counterbalance", json, fixed = TRUE) &&
          grepl("\"spec_version\": \"0.3\"", json, fixed = TRUE),
        "browser app: a between design leaves out the hidden, ticked Counterbalanced lists")
})

# The power curve extends its grid until the grid brackets the 0.80 target, by the rule
# tests/test_app_server.R sets out, except that the browser build stops at 2000 subjects where 8 N
# is more, to stay responsive. The grid used to stop at twice N subjects, which for the default
# design is the target itself, so about half of all runs were refused.
extensions <- list(
  list(what = "a power below the target throughout", power = "below", n = 64,
       grid = c(32, 48, 64, 96, 128, 192, 288, 432, 512),
       says = "does not reach a power of 0.8"),
  list(what = "a power above the target throughout", power = "above", n = 64,
       grid = c(4, 6, 9, 14, 21, 32, 48, 64, 96, 128),
       says = "does not reach a power of 0.8"),
  list(what = "a power that crosses the target above twice N", power = "crossing", n = 64,
       grid = c(32, 48, 64, 96, 128, 192),
       says = "Target power 0.80: N = 126 subjects"),
  list(what = "a power below the target from N = 300", power = "below", n = 300,
       grid = c(150, 225, 300, 450, 600, 900, 1350, 2000),
       says = "does not reach a power of 0.8"))
for (case in extensions) testServer(app = stage, {
  do.call(session$setInputs, start_inputs)
  session$setInputs(n_subject = case$n)
  calls <- new.env()
  restore <- stub_power(stand_ins[[case$power]], calls)
  tryCatch(session$setInputs(run_curve = 1), finally = restore())
  pc <- power_plot()
  text <- shown(output$power_out)
  check(identical(pc$grid, case$grid) && grepl(case$says, text, fixed = TRUE),
        sprintf("browser app: %s gives the grid %s and the line '%s'", case$what,
                paste(case$grid, collapse = ", "), case$says))
  check(identical(calls$seed, replicate_seeds(2024, 10)[seq_along(calls$n)]) &&
          identical(pc$seeds, expected_seeds(2024, case$n, pc$grid)),
        sprintf("browser app: with %s, each point takes the next seed in the order computed",
                case$what))
})

# The analytic sample size is shown for the plain default design, and not for a design with random
# effects or continuous predictors, since either adds variance that the `sd` of power.t.test()
# leaves out.
testServer(app = stage, {
  do.call(session$setInputs, start_inputs)
  by_subject <- current_spec()
  by_subject$random <- list(subject = list(intercept_sd = 5))
  with_age <- current_spec()
  with_age$predictors <- list(list(name = "age", varies_by = "subject", mean = 0, sd = 1))
  with_age$fixed$coefficients$age <- 1
  designs <- list("the default design" = NULL, "a by-subject intercept" = by_subject,
                  "a continuous predictor" = with_age)
  for (i in seq_along(designs)) {
    if (!is.null(designs[[i]]))
      session$setInputs(use_pasted = TRUE, pasted = spec_json(designs[[i]]))
    calls <- new.env()
    restore <- stub_power(stand_ins$crossing, calls)
    tryCatch(session$setInputs(run_curve = i), finally = restore())
    text <- shown(output$power_out)
    shows <- grepl(sprintf("Analytic reference: N = %.1f subjects", analytic_default), text,
                   fixed = TRUE)
    check(grepl("Target power 0.80: N = ", text, fixed = TRUE) && shows == (i == 1L),
          sprintf("browser app: a curve for %s %s the analytic sample size", names(designs)[i],
                  if (i == 1L) "shows" else "solves without"))
  }
})

# With the engine itself, at the 200 simulations a point the browser build runs, the default design
# solves at each of six seeds, and the mean of the six is within 5% of the analytic 127.5 subjects.
solves <- new.env()
for (seed in c(2024, 1:5)) testServer(app = stage, {
  do.call(session$setInputs, start_inputs)
  session$setInputs(seed = seed, n_sims = 300, run_curve = 1)
  pc <- power_plot()
  solves[[as.character(seed)]] <- if (is.null(pc$solved)) NA_real_ else pc$solved$n
  cat(sprintf("  seed %d: %s\n", seed, gsub("\n", " | ", shown(output$power_out))))
})
solved_n <- unlist(mget(as.character(c(2024, 1:5)), envir = solves))
check(!anyNA(solved_n) && abs(mean(solved_n) / analytic_default - 1) < 0.05,
      sprintf("browser app: the default design solves at seeds 2024 and 1 to 5 (N = %s), to a mean within 5%% of %.1f",
              paste(solved_n, collapse = ", "), analytic_default))

cat(if (ok) "LITE TESTSERVER OK\n" else "LITE TESTSERVER FAILED\n")
quit(status = if (ok) 0 else 1)
