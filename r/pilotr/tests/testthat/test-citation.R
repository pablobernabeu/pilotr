# The R and Python packages are released separately, so their version numbers can
# differ (R 0.3.1 reached CRAN while PyPI stayed at 0.3.0). The citation built from
# DESCRIPTION's version must therefore name the R package alone: "R and Python
# package version 0.3.1" cited a Python release that does not exist. Twinned with
# python/tests/test_release_metadata.py, which checks the Python About page.

citation_file <- function() {
  path <- system.file("CITATION", package = "pilotr")
  testthat::skip_if(!nzchar(path), "inst/CITATION not found")
  path
}

test_that("citation('pilotr') names the R package version alone", {
  cit <- utils::readCitationFile(citation_file(), meta = list(Version = "0.3.1"))
  expect_identical(cit$note, "R package version 0.3.1")
  text <- paste(format(cit, style = "textVersion"), collapse = " ")
  expect_match(text, "R package version 0.3.1", fixed = TRUE)
  expect_false(grepl("Python", text, fixed = TRUE))
})

test_that("citation('pilotr') gives the author's ORCID iD", {
  cit <- utils::readCitationFile(citation_file(), meta = list(Version = "0.3.1"))
  author <- cit$author[[1]]
  expect_identical(unname(author$comment["ORCID"]), "0000-0003-1083-2460")
})

test_that("the About article cites the R package version alone", {
  path <- testthat::test_path("..", "..", "vignettes", "articles", "about.Rmd")
  testthat::skip_if_not(file.exists(path), "vignettes/articles/about.Rmd not available")
  text <- readLines(path, encoding = "UTF-8")
  expect_false(any(grepl("R and Python package version", text, fixed = TRUE)))
  expect_identical(sum(grepl("(R package version %s)", text, fixed = TRUE)), 1L)
  expect_identical(sum(grepl("note   = {R package version %s}", text, fixed = TRUE)), 1L)
})
