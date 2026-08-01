# Core public-API tests: spec building, simulation, the JSON round-trip, power, and the
# response families. Kept fast (small n_sims, no lme4) so they run well within CRAN limits.

gaussian_between <- function(seed = 2024) {
  build_spec(list(
    name = "t", seed = seed, n_subject = 64, design_kind = "between",
    factor_name = "group", lev1 = "control", lev2 = "treatment",
    intercept = 100, effect = 5, family = "gaussian", resp_name = "", sigma = 10))
}

test_that("build_spec and simulate_design produce a well-formed between-subjects design", {
  spec <- gaussian_between()
  d <- simulate_design(spec)
  expect_s3_class(d, "data.frame")
  expect_equal(nrow(d), 64)
  expect_true(all(c("subject", "group", "score") %in% names(d)))
  expect_null(spec$random)
  expect_setequal(unique(d$group), c("control", "treatment"))
})

test_that("the same specification and seed reproduce identical data", {
  a <- simulate_design(gaussian_between(seed = 7))
  b <- simulate_design(gaussian_between(seed = 7))
  expect_identical(a, b)
})

test_that("spec_json is valid JSON and round-trips through load_spec", {
  spec <- gaussian_between()
  js <- spec_json(spec)
  expect_true(jsonlite::validate(js))
  tmp <- tempfile(fileext = ".json"); on.exit(unlink(tmp))
  writeLines(js, tmp)
  expect_equal(simulate_design(load_spec(tmp))$score, simulate_design(spec)$score)
})

test_that("power_design returns Type S / Type M and a plausible power", {
  r <- power_design(gaussian_between(), n_sims = 100)
  expect_true(all(c("n_sims", "power", "type_s", "type_m", "true_effect", "mean_estimate")
                  %in% names(r)))
  expect_gte(r$power, 0); expect_lte(r$power, 1)
  expect_gt(r$power, 0.2); expect_lt(r$power, 0.8)   # d = 0.5, n = 32/group
  expect_equal(r$true_effect, 5)
})

test_that("default_response_name covers every family", {
  expect_equal(default_response_name("gaussian"), "score")
  expect_equal(default_response_name("lognormal"), "RT")
  expect_equal(default_response_name("shifted_lognormal"), "RT")
  expect_equal(default_response_name("bernoulli"), "accuracy")
  expect_equal(default_response_name("poisson"), "count")
  expect_equal(default_response_name("ordinal"), "rating")
  expect_equal(default_response_name("beta"), "proportion")
})

test_that("each response family simulates a column on the expected scale", {
  fam <- function(family, extra = list()) {
    p <- c(list(name = "f", seed = 1, n_subject = 80, design_kind = "between",
                factor_name = "group", lev1 = "a", lev2 = "b", intercept = 0,
                effect = 0.5, family = family, resp_name = ""), extra)
    simulate_design(build_spec(p))
  }
  expect_true(all(fam("bernoulli")$accuracy %in% c(0, 1)))
  pois <- fam("poisson", list(intercept = 1.5, effect = 0.3))$count
  expect_true(all(pois >= 0) && all(pois == round(pois)))
  ord <- fam("ordinal", list(intercept = 0, effect = 0.8,
                             thresholds = "-2, -0.6, 0.6, 2"))$rating
  expect_true(all(ord >= 1) && all(ord <= 5))
  prop <- fam("beta", list(intercept = 0, effect = 0.8, phi = 8))$proportion
  expect_true(all(prop > 0) && all(prop < 1))
})

test_that("build_spec carries sigma through for the lognormal family", {
  spec <- build_spec(list(name = "ln", seed = 1, design_kind = "between",
                          factor_name = "g", lev1 = "a", lev2 = "b", n_subject = 10,
                          intercept = 0, effect = 0.5, family = "lognormal",
                          resp_name = "", sigma = 1))
  expect_equal(spec$response$sigma, 1)
  expect_equal(spec$response$name, "RT")
  d <- simulate_design(spec)
  expect_equal(nrow(d), 10)
  expect_true(all(d$RT > 0))
})

test_that("per_subject must lie between 1 and the number of items", {
  spec <- build_spec(list(name = "pc", seed = 1, design_kind = "within",
                          include_items = TRUE, n_subject = 2, n_item = 3,
                          factor_name = "cond", lev1 = "a", lev2 = "b",
                          intercept = 6, effect = 0.05, subj_int_sd = 0.1,
                          subj_slope_sd = 0, item_int_sd = 0.1, item_slope_sd = 0,
                          family = "gaussian", resp_name = "", sigma = 0.3))
  spec$units$item$per_subject <- 5
  expect_error(simulate_design(spec), "cannot exceed the number of items")
  spec$units$item$per_subject <- 0
  expect_error(simulate_design(spec), "at least 1")
})

# Until 0.3 power_mixed() refused any design without an item unit, and any design with more than one
# within-unit factor, because it fitted a formula written into the source rather than the one the
# specification implies. Both restrictions are gone, so what is checked now is that such a design
# runs and reports against the right effects.
test_that("power_mixed accepts a design with no item unit", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  spec <- build_spec(list(name = "w", seed = 1, design_kind = "within",
                          include_items = FALSE, n_subject = 40,
                          factor_name = "cond", lev1 = "a", lev2 = "b",
                          intercept = 6, effect = 0.4, subj_int_sd = 0.12,
                          subj_slope_sd = 0, subj_corr = 0,
                          family = "gaussian", resp_name = "", sigma = 0.3))
  out <- power_mixed(spec, n_sims = 4)
  expect_s3_class(out, "pilotr_power")
  expect_identical(names(out$power), "effect")
  expect_true(out$n_returned > 0L)
})

test_that("power_mixed accepts more than one within-unit factor", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  spec <- build_spec(list(name = "two", seed = 1, design_kind = "within",
                          include_items = TRUE, n_subject = 12, n_item = 8,
                          factor_name = "cond", lev1 = "a", lev2 = "b",
                          intercept = 6, effect = 0.05, subj_int_sd = 0.12,
                          subj_slope_sd = 0, item_int_sd = 0.08, item_slope_sd = 0,
                          family = "gaussian", resp_name = "", sigma = 0.3))
  spec$factors[[2]] <- list(name = "block", levels = c("x", "y"),
                            contrasts = list(blk = c(-0.5, 0.5)),
                            vary_within = c("subject"))
  spec$fixed$coefficients$blk <- 0.03
  out <- power_mixed(spec, n_sims = 4)
  expect_setequal(names(out$power), c("effect", "blk"))
})

test_that("every reported rate carries its Monte Carlo error", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  spec <- build_spec(list(name = "mc", seed = 1, design_kind = "within",
                          include_items = TRUE, n_subject = 10, n_item = 8,
                          factor_name = "cond", lev1 = "a", lev2 = "b",
                          intercept = 6, effect = 0.05, subj_int_sd = 0.12,
                          subj_slope_sd = 0, item_int_sd = 0.08, item_slope_sd = 0,
                          family = "gaussian", resp_name = "", sigma = 0.3))
  out <- power_mixed(spec, n_sims = 6)
  expect_equal(unname(out$power_mcse),
               unname(sqrt(out$power * (1 - out$power) / out$n_returned)), tolerance = 1e-12)
  expect_true(out$power_lo <= out$power && out$power <= out$power_hi)
  # The plain standard error is zero at the boundary, where the Wilson interval still has width.
  expect_identical(pilotr:::.wilson(0, 60)[1], 0)
  expect_true(pilotr:::.wilson(0, 60)[2] > 0)
  expect_identical(pilotr:::.wilson(1, 60)[2], 1)
  expect_true(pilotr:::.wilson(1, 60)[1] < 1)
})

test_that("the fit counts distinguish singular fits from clean convergence", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  # A maximal crossed structure at a small sample size is singular in most replicates, which is the
  # case the old n_converged reported as fully converged.
  spec <- build_spec(list(name = "sing", seed = 1, design_kind = "within",
                          include_items = TRUE, n_subject = 8, n_item = 6,
                          factor_name = "cond", lev1 = "a", lev2 = "b",
                          intercept = 6, effect = 0.05, subj_int_sd = 0.12,
                          subj_slope_sd = 0.04, subj_corr = 0.2,
                          item_int_sd = 0.08, item_slope_sd = 0.02, item_corr = -0.1,
                          family = "shifted_lognormal", resp_name = "", sigma = 0.3, shift = 200))
  out <- power_mixed(spec, n_sims = 6)
  expect_identical(out$n_attempted, 6L)
  expect_true(out$n_converged <= out$n_returned)
  expect_true(out$n_singular > 0L)
})

test_that("an interaction random slope reaches the linear predictor", {
  # The confound-free test: the slope KEY stays present in both runs and only its standard
  # deviation changes, so the number of draws consumed is identical and any difference in the data
  # is attributable to the slope itself. Before 0.3 both runs were bit-identical.
  spec <- load_spec(pilotr_example("reading_time_continuous"))
  spec$spec_version <- "0.3"
  spec$response$round <- NULL
  spec$random$subject$correlations <- NULL
  spec$random$subject$slopes[["SyntaxPC:age"]] <- 0
  zero <- simulate_design(spec)
  spec$random$subject$slopes[["SyntaxPC:age"]] <- 0.4
  nonzero <- simulate_design(spec)
  expect_false(identical(zero$reading_time_per_word, nonzero$reading_time_per_word))
})

test_that("a mistyped coefficient key is refused rather than silently zeroed", {
  spec <- load_spec(pilotr_example("crossed_mixed_rt"))
  spec$fixed$coefficients <- list(cnod = 0.05)
  expect_error(validate_spec(spec), "neither a contrast column nor a predictor")
})

test_that("a non-positive-definite random-effect covariance is an error", {
  spec <- load_spec(pilotr_example("crossed_mixed_rt"))
  spec$random$subject$slopes[["extra"]] <- 0.1
  spec$fixed$coefficients[["extra"]] <- 0
  spec$predictors <- list(list(name = "extra", varies_by = "item", mean = 0, sd = 1))
  spec$random$subject$correlations <- list(`intercept,cond` = 0.95,
                                           `intercept,extra` = 0.95,
                                           `cond,extra` = -0.95)
  expect_error(simulate_design(spec), "not positive definite")
})

test_that("replicate seeds are distinct and reproducible", {
  s <- replicate_seeds(90210, 500)
  expect_length(unique(s), 500)
  expect_identical(s, replicate_seeds(90210, 500))
  expect_true(all(s >= 1))
  # Not the old consecutive rule.
  expect_false(identical(s[1:3], c(90210, 90211, 90212)))
})

test_that("spec_json round-trips a coefficient exactly", {
  spec <- load_spec(pilotr_example("crossed_mixed_rt"))
  spec$fixed$coefficients$cond <- 1 / 3
  f <- tempfile(fileext = ".json")
  writeLines(spec_json(spec), f)
  back <- load_spec(f)
  expect_identical(as.numeric(back$fixed$coefficients$cond), 1 / 3)
  # And it stays readable: a value typed as 0.3 is not written as 0.29999999999999999.
  spec$fixed$coefficients$cond <- 0.3
  expect_true(grepl('"cond": 0.3', spec_json(spec), fixed = TRUE))
})

test_that("brms_bridge maps the beta family to brms's Beta()", {
  spec <- build_spec(list(name = "b", seed = 1, design_kind = "between",
                          factor_name = "g", lev1 = "a", lev2 = "b", n_subject = 10,
                          intercept = 0, effect = 0.8, family = "beta",
                          resp_name = "", phi = 8))
  out <- capture.output(bridge <- brms_bridge(spec))
  expect_equal(bridge$family, "Beta()")
})

test_that("a crossed within-design has by-subject and by-item random slopes", {
  spec <- build_spec(list(
    name = "c", seed = 1, n_subject = 20, n_item = 12, design_kind = "within",
    include_items = TRUE, factor_name = "cond", lev1 = "x", lev2 = "y",
    intercept = 6, effect = 0.05, subj_int_sd = 0.12, subj_slope_sd = 0.04, subj_corr = 0.2,
    item_int_sd = 0.08, item_slope_sd = 0.02, item_corr = -0.1,
    family = "shifted_lognormal", resp_name = "", sigma = 0.3, shift = 200))
  d <- simulate_design(spec)
  expect_equal(nrow(d), 20 * 12 * 2)
  expect_false(is.null(spec$random$subject$slopes))
  expect_false(is.null(spec$random$item$slopes))
  expect_true(all(d$RT > 200))
})

test_that("generate_r_script emits a self-contained, runnable script", {
  rs <- generate_r_script(gaussian_between())
  expect_true(grepl("library(pilotr)", rs, fixed = TRUE))
  expect_true(grepl("simulate_design", rs, fixed = TRUE))
  expect_true(grepl("spec <-", rs, fixed = TRUE))
})

test_that("pilotr_example lists the bundled specs and resolves each to a file", {
  names <- pilotr_example()
  expect_type(names, "character")
  expect_true(length(names) >= 1)
  expect_true("between_2group_gaussian" %in% names)
  for (nm in names) {
    path <- pilotr_example(nm)
    expect_true(file.exists(path))
    # Every shipped example loads and simulates without error.
    spec <- load_spec(path)
    d <- simulate_design(spec)
    expect_s3_class(d, "data.frame")
    expect_gt(nrow(d), 0L)
  }
  # The .json extension is optional.
  expect_identical(
    pilotr_example("between_2group_gaussian"),
    pilotr_example("between_2group_gaussian.json")
  )
})

test_that("pilotr_example rejects unknown or malformed names", {
  expect_error(pilotr_example("no_such_example"), "Unknown example")
  expect_error(pilotr_example(c("a", "b")), "single example name")
})
