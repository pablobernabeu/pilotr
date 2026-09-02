# The power and precision vignettes read three cached results from disk, because the
# mixed-model fits behind them cost far more than a package check may spend. A cache
# outlives the code that wrote it: when a replicate loop changes what it returns, the
# vignette goes on printing the old columns and calling them the output of the code it
# shows. Compare the cached column names against a live, tiny run of the same function,
# so that a change to a return value fails here until the caches are regenerated with
# tools/regenerate-vignette-caches.R.

vignette_cache <- function(file) {
  path <- testthat::test_path("..", "..", "vignettes", file)
  testthat::skip_if_not(file.exists(path), "vignette sources not available")
  path
}

cache_spec <- function(effect) {
  build_spec(list(
    name = "priming", seed = 1, design_kind = "within", include_items = TRUE,
    n_subject = 8, n_item = 6,
    factor_name = "condition", lev1 = "related", lev2 = "unrelated",
    intercept = 6, effect = effect,
    subj_int_sd = 0.12, subj_slope_sd = 0.04, subj_corr = 0.2,
    item_int_sd = 0.08, item_slope_sd = 0.02, item_corr = -0.1,
    family = "shifted_lognormal", resp_name = "RT", sigma = 0.3, shift = 200))
}

test_that("the vignette caches carry the columns the current code returns", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  mixed <- vignette_cache("power-mixed-cache.rds")
  curve <- vignette_cache("power-curve-cache.csv")
  precision <- vignette_cache("precision-curve-cache.csv")
  spec <- cache_spec(0.06)

  pm <- suppressWarnings(power_mixed(spec, n_sims = 2))
  cached_pm <- readRDS(mixed)
  expect_s3_class(cached_pm, "pilotr_power")
  expect_identical(names(cached_pm), names(pm))

  pc <- suppressWarnings(power_curve_mixed(spec, subject_ns = 8, n_sims = 2))
  expect_setequal(names(read.csv(curve)), names(pc))

  prc <- suppressWarnings(precision_curve(cache_spec(0.05), focal = c(effect = 0.05),
                                          subject_ns = 8, rope = 0.02, n_sims = 2))
  expect_setequal(names(read.csv(precision)), names(prc))
})
