# The example specifications under the repository's spec/examples/ are the
# canonical copies; both packages carry a mirror so an installed copy can reach
# them. Nothing enforced the mirror, and a load-and-simulate test cannot: a stale
# packaged copy still loads and simulates perfectly well, it simply describes a
# different design from the one the repository documents. Compare the bytes.
# Twinned with python/tests/test_examples.py. Repository-level checks skip
# gracefully when the package is checked in isolation.

repo_examples_dir <- function() {
  path <- testthat::test_path("..", "..", "..", "..", "spec", "examples")
  testthat::skip_if_not(dir.exists(path), "repository specifications not available")
  path
}

read_raw <- function(path) readBin(path, "raw", file.size(path))

test_that("pilotr_example lists and resolves every packaged specification", {
  names <- pilotr_example()
  expect_type(names, "character")
  expect_true("between_2group_gaussian" %in% names)
  for (nm in names) expect_true(file.exists(pilotr_example(nm)))
  # The .json extension is optional.
  expect_identical(pilotr_example("between_2group_gaussian"),
                   pilotr_example("between_2group_gaussian.json"))
})

test_that("pilotr_example rejects an unknown name", {
  expect_error(pilotr_example("no_such_example"), "Unknown example")
})

test_that("packaged examples are byte-identical to the repository copies", {
  canonical <- repo_examples_dir()
  expected <- sort(basename(list.files(canonical, pattern = "\\.json$", full.names = TRUE)))
  expect_gt(length(expected), 0)
  packaged <- sort(paste0(pilotr_example(), ".json"))
  expect_identical(packaged, expected)
  for (nm in expected) {
    expect_identical(read_raw(pilotr_example(nm)),
                     read_raw(file.path(canonical, nm)),
                     info = nm)
  }
})
