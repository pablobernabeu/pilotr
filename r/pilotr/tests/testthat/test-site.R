# The documents and site files around the package. Both apps' guides said that crossed
# mixed-effects power ran through lme4 in Python as well as in R, where Python fits with
# statsmodels. The vignette on the apps named the browser build's tabs alone. The contributing
# guide named a CSV comparison with a tolerance as the check that keeps the twins identical, where
# CI gates on the harness in tools/parity. The note on advanced designs predated 0.3. The R site
# linked three favicon files that pkgdown/favicon did not hold, the browser build's page was
# titled "Shiny App", and the hub's theme toggle stopped working wherever the browser blocks
# site storage. The checks on files outside the installed package skip when the repository is
# not there. Twinned with python/tests/test_site.py.

site_file <- function(...) {
  path <- testthat::test_path("..", "..", ...)
  testthat::skip_if_not(file.exists(path), paste(file.path(...), "not available"))
  path
}

repo_path <- function(...) site_file("..", "..", ...)

# One line of text, so that a sentence wrapped anywhere is still found.
one_line <- function(path) {
  gsub("\\s+", " ", paste(readLines(path, encoding = "UTF-8", warn = FALSE), collapse = " "))
}

# The text an app shows, approximately: its string literals joined in the order they appear, so
# that a sentence split around tags$code() reads as one.
app_strings <- function(path) {
  pd <- utils::getParseData(parse(path, keep.source = TRUE, encoding = "UTF-8"))
  s <- pd[pd$token == "STR_CONST", ]
  s <- s[order(s$line1, s$col1), ]
  # getParseData() abbreviates a long literal, such as the browser build's data URI, in its text
  # column. getParseText() gives the literal in full.
  literals <- utils::getParseText(pd, s$id)
  paste(vapply(literals, function(x) eval(parse(text = x)), character(1)), collapse = "")
}

# The JavaScript of a page with every try block removed, so that what remains is the code that
# runs unguarded.
unguarded_script <- function(path) {
  html <- paste(readLines(path, encoding = "UTF-8", warn = FALSE), collapse = "\n")
  js <- paste(regmatches(html, gregexpr("(?s)<script[^>]*>.*?</script>", html, perl = TRUE))[[1]],
              collapse = "\n")
  js <- gsub("(^|\\s)//[^\n]*", "\\1", js, perl = TRUE)  # line comments
  repeat {
    m <- regexpr("\\btry\\s*\\{", js, perl = TRUE)
    if (m < 0) break
    chars <- strsplit(js, "")[[1]]
    depth <- 0L
    for (i in seq(m + attr(m, "match.length") - 1L, length(chars))) {
      if (chars[i] == "{") depth <- depth + 1L
      if (chars[i] == "}") depth <- depth - 1L
      if (depth == 0L) break
    }
    js <- paste0(substr(js, 1L, m - 1L), substr(js, i + 1L, nchar(js)))
  }
  js
}

installed_app <- function() {
  path <- system.file("app", "app.R", package = "pilotr")
  testthat::skip_if(!nzchar(path), "inst/app/app.R not found")
  path
}

backends_named <- function(path) {
  text <- app_strings(path)
  !grepl("(via lme4) in R and in Python", text, fixed = TRUE) &&
    grepl("in R (lme4, the reference) and in Python (statsmodels)", text, fixed = TRUE)
}

test_that("the installed app's guide names the backend each language fits mixed models with", {
  expect_true(backends_named(installed_app()))
})

test_that("the browser build's guide and the vignette name the same backends", {
  expect_true(backends_named(repo_path("app-lite", "app.R")))
  vignette <- one_line(site_file("vignettes", "the-no-code-app.Rmd"))
  expect_true(grepl("in Python, through `statsmodels`", vignette, fixed = TRUE))
})

test_that("the vignette on the apps names the tabs of both", {
  tabs <- function(path, fun) {
    src <- paste(readLines(path, encoding = "UTF-8", warn = FALSE), collapse = "\n")
    pattern <- paste0(fun, "\\(\\s*\"([^\"]+)\"")
    sub(pattern, "\\1", regmatches(src, gregexpr(pattern, src))[[1]])
  }
  installed <- tabs(installed_app(), "tabPanel")
  lite <- tabs(repo_path("app-lite", "app.R"), "nav_panel")
  expect_true("Reproducible R script" %in% installed)
  expect_true("R script" %in% lite)
  vignette <- one_line(site_file("vignettes", "the-no-code-app.Rmd"))
  for (tab in union(installed, lite))
    expect_true(grepl(paste0("'", tab, "'"), vignette, fixed = TRUE), info = tab)
})

test_that("the browser build's guide points to the help its paste box carries", {
  # It said that "the example buttons there illustrate the format", and the app has no such
  # buttons. The paste box links to the example specifications and to the format.
  text <- app_strings(repo_path("app-lite", "app.R"))
  expect_false(grepl("example buttons", text, fixed = TRUE))
  expect_true(grepl("Ready-to-run example specifications, and the format itself, are linked",
                    text, fixed = TRUE))
  expect_true(grepl("Ready-to-run example specs are in", text, fixed = TRUE))
})

test_that("the installed app's paste box has a label, as the browser build's has", {
  # textAreaInput() with a NULL label renders an empty <label>, which names the control for
  # nobody, so a screen reader had only the placeholder to read out.
  skip_if_not_installed("shiny")
  label_of <- function(path, id) {
    pd <- utils::getParseData(parse(path, keep.source = TRUE, encoding = "UTF-8"))
    heads <- pd$parent[pd$token == "SYMBOL_FUNCTION_CALL" & pd$text == "textAreaInput"]
    calls <- lapply(utils::getParseText(pd, pd$parent[match(heads, pd$id)]), str2lang)
    calls <- Filter(function(cl) identical(match.call(shiny::textAreaInput, cl)$inputId, id),
                    calls)
    expect_length(calls, 1L)
    match.call(shiny::textAreaInput, calls[[1]])$label
  }
  installed <- label_of(installed_app(), "spec_json_in")
  expect_true(is.character(installed) && nzchar(installed))
  expect_identical(installed, label_of(repo_path("app-lite", "app.R"), "pasted"))
})

test_that("the contributing guides name the parity gates CI runs", {
  gates <- c("tools/parity/run_r.R", "tools/parity/run_py.py", "tools/parity/compare.py",
             "tools/parity/validate_cross.py", "tools/parity/solve_cross.py")
  # The repository's guide and the copy pkgdown renders for the R site.
  for (path in c(repo_path(".github", "CONTRIBUTING.md"),
                 site_file(".github", "CONTRIBUTING.md"))) {
    text <- one_line(path)
    for (gate in gates) expect_true(grepl(gate, text, fixed = TRUE), info = paste(path, gate))
    expect_false(grepl("exercises that contract and must stay green", text, fixed = TRUE),
                 info = path)
    expect_false(grepl("# confirm cross-language agreement", text, fixed = TRUE), info = path)
    expect_true(grepl("maths library", text, fixed = TRUE), info = path)
  }
})

test_that("the R site's copy of the contributing guide is the original below its note", {
  # The copy's note says the two are kept in step by hand, and the copy had fallen behind by two
  # sections.
  original <- readLines(repo_path(".github", "CONTRIBUTING.md"), encoding = "UTF-8", warn = FALSE)
  copy <- readLines(site_file(".github", "CONTRIBUTING.md"), encoding = "UTF-8", warn = FALSE)
  expect_identical(copy[seq(match(original[1], copy), length(copy))], original)
})

test_that("the note on advanced designs describes specification 0.3", {
  text <- one_line(repo_path("docs", "mixed_models_and_design_analysis.md"))
  for (term in c("\"observation\"", "\"uniform\"", "reliability", "exgaussian"))
    expect_true(grepl(term, text, fixed = TRUE), info = term)
  expect_false(grepl("5e-15", text, fixed = TRUE))
})

test_that("pkgdown/favicon holds every favicon pkgdown links, and nothing else", {
  # pkgdown 2.1.2 changed the favicons its pages link, to the set RealFaviconGenerator makes now.
  skip_if_not_installed("pkgdown", "2.1.2")
  dir <- site_file("pkgdown", "favicon")
  head <- system.file("BS5", "templates", "head.html", package = "pkgdown")
  skip_if(!nzchar(head), "pkgdown's page head template not found")
  template <- paste(readLines(head, warn = FALSE), collapse = "\n")
  block <- regmatches(template, regexpr("(?s)\\{\\{#has_favicons\\}\\}.*?\\{\\{/has_favicons\\}\\}",
                                        template, perl = TRUE))
  expect_length(block, 1)
  pattern <- "href=\"\\{\\{#site\\}\\}\\{\\{root\\}\\}\\{\\{/site\\}\\}([^\"]+)\""
  linked <- sub(pattern, "\\1", regmatches(block, gregexpr(pattern, block))[[1]])
  expect_true("favicon-96x96.png" %in% linked)
  expect_true(all(file.exists(file.path(dir, linked))),
              info = paste(linked[!file.exists(file.path(dir, linked))], collapse = ", "))
  manifest <- jsonlite::fromJSON(file.path(dir, "site.webmanifest"), simplifyVector = FALSE)
  icons <- vapply(manifest$icons, function(i) sub("^/", "", i$src), character(1))
  expect_true(length(icons) > 0 && all(file.exists(file.path(dir, icons))))
  expect_setequal(list.files(dir), union(linked, icons))
})

test_that("the browser build's page has the app's title and an icon", {
  for (path in c(repo_path(".github", "workflows", "site.yml"),
                 repo_path("app-lite", "build_shinylive.R"))) {
    text <- one_line(path)
    expect_true(grepl("shinylive::export(", text, fixed = TRUE), info = path)
    call <- sub("^.*?shinylive::export\\(", "", text, perl = TRUE)
    expect_true(grepl("template_params", call, fixed = TRUE), info = path)
    expect_true(grepl("title = \"pilotr: no-code app\"", text, fixed = TRUE), info = path)
    expect_true(grepl("include_in_head", call, fixed = TRUE), info = path)
    expect_true(grepl("rel=\"icon\"", text, fixed = TRUE), info = path)
    # The icon the page links is copied beside it.
    expect_true(grepl("favicon.svg", text, fixed = TRUE), info = path)
  }
})

test_that("the hub's pages read and write site storage inside try blocks", {
  for (page in c("landing.html", "404.html")) {
    path <- repo_path(".github", page)
    expect_true(grepl("localStorage", one_line(path), fixed = TRUE), info = page)
    expect_false(grepl("localStorage", unguarded_script(path), fixed = TRUE), info = page)
  }
})
