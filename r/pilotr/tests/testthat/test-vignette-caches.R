# The power and precision vignettes show three chunks with eval = FALSE and read their results from
# files in vignettes/. The hundreds of mixed-model fits behind them take far longer than a package
# check allows. A cached result outlives the code that produced it, though. Up to 0.3.1, the
# vignettes printed results from before 0.3.0, which had changed the replicate seeds and the fit
# accounting. They presented those results as the output of the chunk above them, and lacked
# columns the analyses had returned since. These tests fail when the columns, the seeds, the
# simulated data or the chunks behind the caches change, until tools/vignette-caches.R writes the
# caches afresh.
#
# Both need the package's source tree, where vignettes/ and tools/ sit beside tests/. They skip
# when vignettes/ is absent, as it is in an installed package and under R CMD check, and on CRAN.
# The second also skips without tools/, which the built package leaves out.

vignette_dir <- function() {
  dir <- testthat::test_path("..", "..", "vignettes")
  skip_if_not(dir.exists(dir), "the vignette sources are not available")
  dir
}

test_that("the cached results have every column the analyses now return", {
  skip_on_cran()
  dir <- vignette_dir()
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  # The crossed design of the vignettes, shrunk so that each analysis runs in a second or two. Only
  # the names of what comes back are compared.
  spec <- build_spec(list(
    name = "priming", seed = 1, design_kind = "within", include_items = TRUE,
    n_subject = 8, n_item = 6,
    factor_name = "condition", lev1 = "related", lev2 = "unrelated",
    intercept = 6, effect = 0.05,
    subj_int_sd = 0.12, subj_slope_sd = 0.04, subj_corr = 0.2,
    item_int_sd = 0.08, item_slope_sd = 0.02, item_corr = -0.1,
    family = "shifted_lognormal", resp_name = "RT", sigma = 0.3, shift = 200))
  # The two tables end with the versions that produced them.
  versions <- c("pilotr_version", "lme4_version")

  curve <- utils::read.csv(file.path(dir, "power-curve-cache.csv"))
  expect_identical(names(curve),
                   c(names(power_curve_mixed(spec, subject_ns = 8, n_sims = 2)), versions))

  prc <- utils::read.csv(file.path(dir, "precision-curve-cache.csv"))
  expect_identical(names(prc),
                   c(names(precision_curve(spec, focal = c(effect = 0.05), subject_ns = 8,
                                           rope = 0.02, n_sims = 2)), versions))

  # The power_mixed() result is kept whole, so that the vignette prints what a live run prints.
  cached <- readRDS(file.path(dir, "power-mixed-cache.rds"))
  expect_s3_class(cached, "pilotr_power")
  expect_identical(names(cached), names(power_mixed(spec, n_sims = 2)))
})

test_that("the cached results rest on the seeds, data and code the vignettes now have", {
  # A change to the seed rule or the generative engine changes every cached number while leaving
  # the columns as they were, and so does an edit to a chunk whose result is cached. The
  # fingerprint the cache tool writes covers all three, and is recomputed here by the tool's own
  # function, so the two cannot disagree about what is being compared.
  skip_on_cran()
  dir <- vignette_dir()
  tool_path <- testthat::test_path("..", "..", "tools", "vignette-caches.R")
  skip_if_not(file.exists(tool_path), "tools/vignette-caches.R is not available")
  tool <- new.env(parent = environment())
  sys.source(tool_path, envir = tool)
  current <- tool$cache_fingerprint(dir)
  written <- read.dcf(file.path(dir, "cache-fingerprint.txt"))[1, ]
  stale <- names(current)[!(names(current) %in% names(written)) |
                            written[names(current)] != current]
  expect(length(stale) == 0L, sprintf(
    "The vignette caches no longer match the code (%s differs). Regenerate them from r/pilotr with `Rscript tools/vignette-caches.R`, and commit what it writes.",
    paste(stale, collapse = ", ")))
})
