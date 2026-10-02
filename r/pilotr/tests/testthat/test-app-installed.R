# These tests cover the bundled app as the installed package runs it. pilotr::run_app() loads the
# namespace without attaching the package, and shiny sources app.R into an environment whose
# parent is the global environment, so the app cannot find pilotr's functions on the search
# path. Until this was fixed, the app looked there and then for source files that an installed
# package does not have, so the documented launcher stopped before the app started. The
# source-tree scripts in tests/ never met this case, because from the source tree the app
# sources R/ itself.

test_that("the installed app simulates in a session that has not attached pilotr", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("callr")
  skip_if_not_installed("ggplot2")

  # This process has attached pilotr, so the app runs in a fresh one that only loads the
  # namespace, as run_app() does. Under devtools::test(), the code under test is the source
  # tree, while a fresh process loads whichever pilotr is installed. The fresh process
  # therefore reports the copy it loaded, and the test goes on only when that copy is the one
  # under test.
  under_test <- normalizePath(getNamespaceInfo("pilotr", "path"), winslash = "/")
  res <- callr::r(function(under_test) {
    if (!requireNamespace("pilotr", quietly = TRUE)) return(list(loaded = NA_character_))
    loaded <- normalizePath(getNamespaceInfo("pilotr", "path"), winslash = "/")
    if (!identical(loaded, under_test)) return(list(loaded = loaded))
    # testServer() evaluates its block inside a mask of the server's environment, so results
    # leave through an environment, which the block modifies in place.
    out <- new.env()
    shiny::testServer(system.file("app", package = "pilotr"), {
      # The values the app's controls start at.
      session$setInputs(
        name = "my_design", seed = 2024, design_kind = "between", n_subject = 64,
        include_items = TRUE, n_item = 24, factor_name = "group", lev1 = "control",
        lev2 = "treatment", intercept = 100, effect = 5, family = "gaussian",
        resp_name = "", sigma = 10, shift = 200, spec_json_in = "")
      session$setInputs(simulate = 1)
      out$n <- nrow(data())
      out$columns <- names(data())
      # With the package's functions in use, there are no engine files to source, so Verify's
      # clean session loads the installed package.
      session$setInputs(verify_code = 1)
      out$verify <- output$verify_out
    })
    list(loaded = loaded, n = out$n, columns = out$columns, verify = out$verify,
         attached = "package:pilotr" %in% search())
  }, args = list(under_test = under_test))
  skip_if(!identical(res$loaded, under_test),
          "a fresh R process does not load the copy of pilotr under test")

  expect_identical(res$n, 64L)
  expect_true(all(c("subject", "group", "score") %in% res$columns))
  expect_false(res$attached)
  expect_match(res$verify, "^Reproduces")
})

test_that("every pilotr object the app uses is exported", {
  # The installed app binds pilotr's exports and nothing else. An internal helper, called or
  # passed as an argument, and an internal constant would all still work from the source tree,
  # where the app sources every file in R/. They would fail only once the package is
  # installed. Helpers that only the app needs live in app.R.
  app <- system.file("app", "app.R", package = "pilotr")
  skip_if(!nzchar(app), "inst/app/app.R not found")
  parsed <- utils::getParseData(parse(app, keep.source = TRUE, encoding = "UTF-8"))
  used <- unique(parsed$text[parsed$token %in% c("SYMBOL_FUNCTION_CALL", "SYMBOL")])
  internal <- setdiff(ls(asNamespace("pilotr"), all.names = TRUE),
                      getNamespaceExports("pilotr"))
  expect_identical(intersect(used, internal), character(0))
})
