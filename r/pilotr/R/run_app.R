# Installable launcher for the no-code app. Each user runs it locally on their own machine,
# so concurrency is unbounded across users and heavy power runs use the user's own cores.

#' Launch the pilotr no-code app.
#'
#' Runs the bundled Shiny app locally. When the package is installed, the app calls the
#' package's functions directly; running locally means each user has a private R process
#' (no shared-process blocking) and can parallelise power simulations across their own cores.
#'
#' @param ... passed to [shiny::runApp()] (e.g. `port`, `launch.browser`).
#' @param async If TRUE (default when 'future' and 'promises' are installed), set a
#'   multisession future plan so power runs execute in a background worker and never block
#'   the UI -- the relevant case if you deploy the app as a shared multi-user instance.
#' @export
run_app <- function(..., async = NULL) {
  if (!requireNamespace("shiny", quietly = TRUE)) stop("Please install 'shiny' to run the app.")
  if (is.null(async)) async <- requireNamespace("future", quietly = TRUE) &&
      requireNamespace("promises", quietly = TRUE)
  if (isTRUE(async)) {
    future::plan(future::multisession)
    on.exit(future::plan(future::sequential), add = TRUE)
  }
  app_dir <- system.file("app", package = "pilotr")
  if (!nzchar(app_dir)) {                      # dev fallback when running from source
    cand <- c("inst/app", file.path("toolkit", "r", "pilotr", "inst", "app"))
    app_dir <- cand[file.exists(file.path(cand, "app.R"))][1]
    if (is.na(app_dir)) stop("Could not locate the app directory.")
  }
  shiny::runApp(app_dir, ...)
}
