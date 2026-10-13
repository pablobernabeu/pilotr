# Parallel execution must reproduce the serial stream exactly. Every replicate takes its own
# seed from replicate_seeds(), so the worker count can change only the wall-clock time,
# never a number. Cluster tests are skipped on CRAN to respect its two-core limit and to
# avoid PSOCK flakiness on the check farms.

gaussian_between_par <- function(seed = 11) {
  build_spec(list(
    name = "par", seed = seed, n_subject = 24, design_kind = "between",
    factor_name = "group", lev1 = "control", lev2 = "treatment",
    intercept = 100, effect = 5, family = "gaussian", resp_name = "", sigma = 10))
}

crossed_within_par <- function() {
  build_spec(list(
    name = "parmix", seed = 3, n_subject = 12, n_item = 8, design_kind = "within",
    include_items = TRUE, factor_name = "cond", lev1 = "a", lev2 = "b",
    intercept = 6, effect = 0.05, subj_int_sd = 0.12, subj_slope_sd = 0.04, subj_corr = 0.2,
    item_int_sd = 0.08, item_slope_sd = 0.02, item_corr = -0.1,
    family = "shifted_lognormal", resp_name = "", sigma = 0.3, shift = 200))
}

test_that("workers is validated as a single positive whole number", {
  spec <- gaussian_between_par()
  expect_error(power_design(spec, n_sims = 2, workers = 0), "workers")
  expect_error(power_design(spec, n_sims = 2, workers = 1.5), "workers")
  expect_error(power_design(spec, n_sims = 2, workers = c(2, 4)), "workers")
  expect_error(power_design(spec, n_sims = 2, workers = "two"), "workers")
})

test_that("power_design with two workers is identical to a serial run", {
  skip_on_cran()
  spec <- gaussian_between_par()
  serial <- power_design(spec, n_sims = 8, workers = 1)
  parallel2 <- power_design(spec, n_sims = 8, workers = 2)
  expect_identical(serial, parallel2)
})

test_that("power_mixed and power_curve_mixed with two workers are identical to serial runs", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  spec <- crossed_within_par()
  serial <- power_mixed(spec, n_sims = 4, workers = 1)
  parallel2 <- power_mixed(spec, n_sims = 4, workers = 2)
  expect_identical(serial, parallel2)

  curve_serial <- power_curve_mixed(spec, subject_ns = c(12, 14), n_sims = 3, workers = 1)
  curve_parallel <- power_curve_mixed(spec, subject_ns = c(12, 14), n_sims = 3, workers = 2)
  expect_identical(curve_serial, curve_parallel)
})

test_that("precision_design with two workers is identical to a serial run", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  spec <- crossed_within_par()
  serial <- precision_design(spec, focal = c(effect = 0.05), rope = 0.02,
                             n_sims = 4, workers = 1)
  parallel2 <- precision_design(spec, focal = c(effect = 0.05), rope = 0.02,
                                n_sims = 4, workers = 2)
  expect_identical(serial, parallel2)
})

# Whether PSOCK workers load the copy of pilotr under test. Under devtools::test() that copy is the
# source tree, while a worker loads whichever pilotr is installed, and an older installation lacks
# internals that the replicates call. A worker is asked, since it sees the library paths that
# workers see.
workers_load_this_pilotr <- function() {
  under_test <- normalizePath(getNamespaceInfo("pilotr", "path"), winslash = "/")
  probe <- function() {
    if (!requireNamespace("pilotr", quietly = TRUE)) return(NA_character_)
    normalizePath(getNamespaceInfo("pilotr", "path"), winslash = "/")
  }
  environment(probe) <- globalenv()
  cl <- parallel::makeCluster(1L)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  identical(parallel::clusterCall(cl, probe)[[1L]], under_test)
}
skip_unless_workers_load_this_pilotr <- function() {
  skip_if_not(workers_load_this_pilotr(),
              "PSOCK workers load a pilotr other than the copy under test")
}

# The value of `expr` and the number of PSOCK clusters started while it ran. The count lives in an
# environment written into the tracer, since the tracer runs in parallel's namespace.
count_clusters <- function(expr) {
  made <- new.env()
  made$n <- 0L
  suppressMessages(trace("makeCluster", where = asNamespace("parallel"), print = FALSE,
    tracer = bquote(assign("n", get("n", envir = .(made)) + 1L, envir = .(made)))))
  on.exit(suppressMessages(untrace("makeCluster", where = asNamespace("parallel"))), add = TRUE)
  value <- expr
  list(value = value, clusters = made$n)
}

test_that("a curve over sample size starts one worker pool for the whole sweep", {
  # Every grid point used to start and stop a pool of its own, which on a loaded machine cost
  # tens of seconds a point, while the documentation said that one pool was reused.
  skip_on_cran()
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  skip_unless_workers_load_this_pilotr()
  spec <- crossed_within_par()
  mixed <- count_clusters(power_curve_mixed(spec, subject_ns = c(12, 14), n_sims = 2,
                                            workers = 2))
  expect_identical(mixed$clusters, 1L)
  expect_identical(mixed$value, power_curve_mixed(spec, subject_ns = c(12, 14), n_sims = 2))
  precision <- count_clusters(precision_curve(spec, focal = c(effect = 0.05),
                                              subject_ns = c(12, 14), rope = 0.02, n_sims = 2,
                                              workers = 2))
  expect_identical(precision$clusters, 1L)
  expect_identical(precision$value,
                   precision_curve(spec, focal = c(effect = 0.05), subject_ns = c(12, 14),
                                   rope = 0.02, n_sims = 2))
})

test_that("an analysis given a cluster runs on it and leaves it running", {
  skip_on_cran()
  skip_unless_workers_load_this_pilotr()
  spec <- gaussian_between_par()
  cl <- parallel::makeCluster(2L)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  run <- count_clusters(power_design(spec, n_sims = 8, cl = cl))
  expect_identical(run$clusters, 0L)
  expect_identical(run$value, power_design(spec, n_sims = 8))
  expect_identical(parallel::clusterCall(cl, function() 1L), list(1L, 1L))
  # An effect-size sweep hands the same cluster to every grid point through sweep_spec()'s `...`.
  sweep <- count_clusters(sweep_spec(spec, "fixed$coefficients$effect", c(2, 4), power_design,
                                     n_sims = 8, cl = cl))
  expect_identical(sweep$clusters, 0L)
  expect_identical(sweep$value, sweep_spec(spec, "fixed$coefficients$effect", c(2, 4),
                                           power_design, n_sims = 8))
})

test_that("cl must be a cluster", {
  refusal <- "`cl` must be a cluster made by parallel::makeCluster(), or NULL."
  expect_error(power_design(gaussian_between_par(), n_sims = 2, cl = 2), refusal, fixed = TRUE)
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  spec <- crossed_within_par()
  expect_error(power_mixed(spec, n_sims = 2, cl = "cluster"), refusal, fixed = TRUE)
  expect_error(precision_design(spec, n_sims = 2, cl = list()), refusal, fixed = TRUE)
})
