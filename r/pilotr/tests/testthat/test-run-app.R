# run_app() sets a future plan only when the caller has none, puts back the plan it found, and
# tells the app through the option pilotr.app_async whether to run power in a background R
# process. It used to set a multisession plan with a worker for every core whenever `async` was
# TRUE and to end by setting a sequential plan, which discarded a plan the caller had set. The app
# ran power in a background process whenever future and promises were installed, whatever `async`
# said.

skip_unless_app_async <- function() {
  skip_on_cran()
  for (pkg in c("shiny", "ggplot2", "future", "promises")) skip_if_not_installed(pkg)
}

# The option and the future plan in force while run_app() calls shiny::runApp(), read by a tracer
# that then stops runApp() before it starts anything. The tracer runs in shiny's namespace, so what
# it reads is written into an environment.
launch_state <- function(...) {
  seen <- new.env()
  suppressMessages(trace("runApp", where = asNamespace("shiny"), print = FALSE,
    tracer = bquote({
      assign("async", getOption("pilotr.app_async"), envir = .(seen))
      assign("plan", class(future::plan()), envir = .(seen))
      assign("workers", future::nbrOfWorkers(), envir = .(seen))
      stop("runApp() stopped by the test")
    })))
  on.exit(suppressMessages(untrace("runApp", where = asNamespace("shiny"))), add = TRUE)
  expect_error(run_app(...), "stopped by the test")
  as.list(seen)
}

# testServer() sources app.R, whose library() calls attach shiny and ggplot2. Whatever they attach
# is detached again, so that the test files after this one see the search path they would have
# seen without it.
with_app <- function(code) {
  before <- search()
  on.exit(for (name in setdiff(search(), before)) detach(name, character.only = TRUE),
          add = TRUE)
  code
}

test_that("run_app() keeps a future plan the caller has set", {
  skip_unless_app_async()
  oplan <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(oplan), add = TRUE)
  inside <- launch_state(async = TRUE)
  expect_true("multisession" %in% inside$plan)
  expect_equal(inside$workers, 2)
  expect_true(inherits(future::plan(), "multisession"))
  expect_equal(future::nbrOfWorkers(), 2)
})

test_that("run_app() sets a two-worker plan over a sequential one and restores both settings", {
  skip_unless_app_async()
  oplan <- future::plan(future::sequential)
  on.exit(future::plan(oplan), add = TRUE)
  oopt <- options(pilotr.app_async = "the caller's")
  on.exit(options(oopt), add = TRUE)
  inside <- launch_state(async = TRUE)
  expect_true("multisession" %in% inside$plan)
  expect_equal(inside$workers, 2)
  expect_identical(inside$async, TRUE)
  expect_true(inherits(future::plan(), "sequential"))
  expect_identical(getOption("pilotr.app_async"), "the caller's")
})

test_that("run_app(async = FALSE) sets no plan and has the app run power in its own process", {
  skip_unless_app_async()
  oplan <- future::plan(future::sequential)
  on.exit(future::plan(oplan), add = TRUE)
  oopt <- options(pilotr.app_async = "the caller's")
  on.exit(options(oopt), add = TRUE)
  inside <- launch_state(async = FALSE)
  expect_identical(inside$async, FALSE)
  expect_true("sequential" %in% inside$plan)
  expect_identical(getOption("pilotr.app_async"), "the caller's")
})

test_that("the app runs power in a background process only when run_app() asks it to", {
  skip_unless_app_async()
  app <- system.file("app", package = "pilotr")
  skip_if(!nzchar(app), "inst/app not found")
  seen <- new.env()
  oopt <- options(pilotr.app_async = NULL)
  on.exit(options(oopt), add = TRUE)
  # The block testServer() evaluates sees the session's own objects, while `.async_ok` is set once
  # for the app, in the environment enclosing the server function.
  with_app(for (async in list(NULL, FALSE, TRUE)) {
    options(pilotr.app_async = async)
    label <- if (is.null(async)) "unset" else format(async)
    shiny::testServer(app, assign(label, get(".async_ok", envir = session$env), envir = seen))
  })
  # Without the option, as under shiny::runApp() on the app's directory, power runs in the app's
  # own process.
  expect_false(seen[["unset"]])
  expect_false(seen[["FALSE"]])
  expect_true(seen[["TRUE"]])
})

test_that("with async, the app's power curve runs in a future and gives the same curve", {
  # The curve's grid is extended point by point, so the whole extension runs in one future, which
  # returns the curve for the app to solve. Under a sequential plan the future runs in this process,
  # so the curve is the one computed without async.
  skip_unless_app_async()
  app <- system.file("app", package = "pilotr")
  skip_if(!nzchar(app), "inst/app not found")
  oplan <- future::plan(future::sequential)
  on.exit(future::plan(oplan), add = TRUE)
  oopt <- options(pilotr.app_async = NULL)
  on.exit(options(oopt), add = TRUE)
  out <- new.env()
  with_app(for (async in c(FALSE, TRUE)) {
    options(pilotr.app_async = async)
    shiny::testServer(app, {
      session$setInputs(
        name = "t", seed = 2024, design_kind = "between", n_subject = 32, include_items = FALSE,
        n_item = 24, factor_name = "group", lev1 = "control", lev2 = "treatment",
        intercept = 100, effect = 10, family = "gaussian", resp_name = "", sigma = 10,
        n_sims = 100, spec_json_in = "")
      session$setInputs(run_curve = 1)
      # The future's promise settles on the event loop, which each elapse() runs once, and the
      # outputs are rendered again once it has.
      for (i in seq_len(600)) {
        if (!is.null(power_curve_data())) break
        Sys.sleep(0.05)
        session$elapse(50)
      }
      session$flushReact()
      assign(format(async), list(curve = power_curve_data(), text = output$power_out),
             envir = out)
    })
  })
  expect_false(is.null(out[["TRUE"]]$curve))
  expect_match(out[["TRUE"]]$text, "Target power 0.80: N = ", fixed = TRUE)
  expect_identical(out[["TRUE"]], out[["FALSE"]])
})
