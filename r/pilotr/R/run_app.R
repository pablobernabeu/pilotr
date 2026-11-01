# Installable launcher for the no-code app. Each user runs it locally on their own machine, so
# concurrency is unbounded across users. Heavy power runs belong in the console, where
# power_design() and power_mixed() take a `workers` argument.

#' Launch the pilotr no-code app
#'
#' Runs the bundled Shiny app locally. When the package is installed, the app calls the
#' package's functions directly. Running locally gives each user a private R process, so no one
#' else's work waits on theirs.
#'
#' The app needs the 'shiny' and 'ggplot2' packages, which pilotr only suggests, and its Verify
#' button needs the 'callr' package.
#'
#' Each power analysis and power curve in the app runs in one R process, and the session that
#' asked for it waits for the result. Its controls and outputs do not update until the run ends.
#' For a heavy run, download the specification and call [power_design()] or [power_mixed()] from
#' the console, where `workers` spreads the replicates over your local cores.
#'
#' @param ... passed to [shiny::runApp()] (e.g. `port`,
#'   `launch.browser`).
#' @param async If TRUE (the default when 'future' and 'promises' are installed), the app runs
#'   each power analysis and power curve in a background R process. That leaves the app's own R
#'   process free for other sessions, such as a second browser tab, while the session that
#'   started the run waits for it. A parallel future plan you have set is used as it stands. Over a
#'   sequential plan, the default, `run_app()` sets a multisession plan of two workers while the
#'   app runs. Either way, the plan in force before the call is in force again once the app stops.
#'   With FALSE, power runs in the app's own R process and no plan is set.
#' @return No return value; called for its side effect of launching the Shiny application,
#'   which blocks the R session until the app is closed.
#' @examples
#' \dontrun{
#' run_app()
#' }
#' @export
run_app <- function(..., async = NULL) {
  if (!requireNamespace("shiny", quietly = TRUE)) stop("Please install 'shiny' to run the app.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install 'ggplot2' to run the app.")
  async_installed <- requireNamespace("future", quietly = TRUE) &&
    requireNamespace("promises", quietly = TRUE)
  async <- if (is.null(async)) async_installed else isTRUE(async)
  if (async && !async_installed)
    stop("Please install 'future' and 'promises' to run power in the background, ",
         "or set async = FALSE.")
  # The app reads this option when shiny sources app.R, which happens inside runApp() below, and
  # runs power in a background process only when it is TRUE.
  old_options <- options(pilotr.app_async = async)
  on.exit(options(old_options), add = TRUE)
  # A plan the caller has set is theirs, and is kept. Over a sequential one, two workers are
  # enough: a session runs one power analysis at a time, and each further worker is an idle R
  # process. future's vignette for package developers ("Avoid changing the future backend") asks
  # a function that sets a plan to undo it on exit, on error as well.
  if (async && inherits(future::plan(), "sequential")) {
    old_plan <- future::plan(future::multisession, workers = 2)
    on.exit(future::plan(old_plan), add = TRUE)
  }
  app_dir <- system.file("app", package = "pilotr")
  if (!nzchar(app_dir)) {                      # dev fallback when running from source
    cand <- c("inst/app", file.path("r", "pilotr", "inst", "app"))
    app_dir <- cand[file.exists(file.path(cand, "app.R"))][1]
    if (is.na(app_dir)) stop("Could not locate the app directory.")
  }
  shiny::runApp(app_dir, ...)
}
