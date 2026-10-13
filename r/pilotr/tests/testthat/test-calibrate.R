# calibrate_response() holds the intercept under tune = "all", since the intercept is a location.
# Scaling it with the rest raised every lognormal response to a power, so the documented example's
# median reaction time of 600 ms became 16.6 hours. A poisson design was also refused a target it
# could reach. denominator = "random" calibrates the random components plus the residual, the
# variance Westfall, Kenny and Judd (2014) divide by to define Cohen's d. The figures in the
# comments are what the function returned before these changes.

example_spec <- function(name) load_spec(pilotr_example(name))

# The specification with every coefficient, random-effect standard deviation and residual scale
# multiplied by k, and the intercept left alone, which is what tune = "all" now does.
scaled_by <- function(spec, k) {
  for (nm in names(spec$fixed$coefficients))
    spec$fixed$coefficients[[nm]] <- spec$fixed$coefficients[[nm]] * k
  for (g in names(spec$random)) {
    spec$random[[g]]$intercept_sd <- spec$random[[g]]$intercept_sd * k
    for (nm in names(spec$random[[g]]$slopes))
      spec$random[[g]]$slopes[[nm]] <- spec$random[[g]]$slopes[[nm]] * k
  }
  for (nm in intersect(c("sigma", "beta"), names(spec$response)))
    spec$response[[nm]] <- spec$response[[nm]] * k
  spec
}

# Cohen's d estimated from simulated data: the difference between the group means over the
# pooled within-group standard deviation.
sample_d <- function(spec) {
  d <- simulate_design(spec)
  g <- d$group == "treatment"
  pooled <- sqrt((stats::var(d$score[g]) + stats::var(d$score[!g])) / 2)
  (mean(d$score[g]) - mean(d$score[!g])) / pooled
}

test_that("tune = \"all\" holds the intercept, so reaction times keep their median", {
  # The intercept went from 6 to 17.93, and the median reaction time from 600 ms to 6e7 ms.
  spec <- example_spec("crossed_mixed_rt")
  cal <- calibrate_response(spec, target_var = 1, tune = "all")
  expect_identical(cal$fixed$intercept, 6)
  expect_equal(response_variance(cal)$total, 1, tolerance = 1e-10)
  # Everything that carries variance is scaled by one factor, the one used before.
  k <- cal$response$sigma / spec$response$sigma
  expect_equal(k, 2.9890666968708324, tolerance = 1e-12)
  expect_equal(cal, scaled_by(spec, k), tolerance = 1e-12)
  rt <- simulate_design(cal)$RT
  expect_gt(median(rt), 500)
  expect_lt(median(rt), 720)
})

test_that("a lognormal outcome keeps its intercept, and its median with it", {
  # The intercept went from -1.2 to -3.51, and the median reading time per word from 0.307 to
  # 0.0316.
  spec <- example_spec("reading_time_continuous")
  cal <- calibrate_response(spec, target_var = 1, tune = "all")
  expect_identical(cal$fixed$intercept, -1.2)
  expect_equal(response_variance(cal)$total, 1, tolerance = 1e-10)
  expect_equal(median(simulate_design(cal)$reading_time_per_word), exp(-1.2), tolerance = 0.25)
})

test_that("for an identity link the intercept is the only part of the output that changes", {
  for (nm in c("between_2group_gaussian", "partial_crossing")) {
    spec <- example_spec(nm)
    held <- calibrate_response(spec, 1, tune = "all")
    scaled <- calibrate_response(spec, 1, tune = "all", scale_intercept = TRUE)
    k <- held$response$sigma / spec$response$sigma
    expect_identical(held$fixed$intercept, spec$fixed$intercept, info = nm)
    expect_equal(scaled$fixed$intercept, spec$fixed$intercept * k, info = nm)
    held$fixed$intercept <- scaled$fixed$intercept
    expect_identical(held, scaled, info = nm)
  }
  # The ex-Gaussian's exponential scale is rescaled with sigma, and its intercept held.
  rt <- example_spec("between_2group_gaussian")
  rt$spec_version <- "0.3"
  rt$fixed <- list(intercept = 500, coefficients = list(grp = 30))
  rt$response <- list(family = "exgaussian", name = "score", sigma = 40, beta = 120)
  cal <- calibrate_response(rt, 1, tune = "all")
  expect_identical(cal$fixed$intercept, 500)
  expect_equal(response_variance(cal)$total, 1, tolerance = 1e-10)
  expect_equal(cal, scaled_by(rt, cal$response$beta / 120), tolerance = 1e-12)
})

test_that("a poisson design calibrates to targets of 1 and 2, and the result simulates", {
  # Target 1 was refused with "even with every effect set to zero the poisson family contributes
  # a residual variance of 1.64493". That is the residual at a linear predictor of zero. At this
  # design's intercept it is 0.221, and the design's own total is 0.267. Target 2 returned an
  # intercept of 11.38, a mean count of 21235 that simulate_design() refuses.
  spec <- example_spec("poisson_counts_between")
  for (target in c(1, 2)) {
    cal <- calibrate_response(spec, target, tune = "all")
    expect_identical(cal$fixed$intercept, spec$fixed$intercept)
    expect_equal(response_variance(cal)$total, target, tolerance = 1e-6)
    counts <- simulate_design(cal)$count
    expect_true(all(is.finite(counts) & counts >= 0))
  }
  # At target 1 the factor is 4.087, which puts the treatment group's linear predictor at 2.43.
  cal <- calibrate_response(spec, 1, tune = "all")
  expect_equal(cal$fixed$coefficients$grp / spec$fixed$coefficients$grp, 4.087, tolerance = 1e-3)
})

test_that("a target the intercept alone already reaches is refused with that residual", {
  spec <- example_spec("poisson_counts_between")
  at_intercept <- trigamma(exp(spec$fixed$intercept))   # about 0.221, at a mean count of 5
  expect_error(
    calibrate_response(spec, 0.2, tune = "all"),
    sprintf(paste0("with every coefficient and random-effect standard deviation set to zero, the ",
                   "intercept alone gives the poisson family a residual variance of %.6g on the ",
                   "latent scale, which already meets or exceeds the target of 0.2; raise ",
                   "`target_var` above %.6g"), at_intercept, at_intercept),
    fixed = TRUE)
})

test_that("scale_intercept = TRUE reproduces the earlier output", {
  # What tune = "all" returned before it held the intercept.
  old <- example_spec("crossed_mixed_rt")
  old$fixed$intercept <- 17.934400181224994
  old$fixed$coefficients$cond <- 0.14945333484354162
  old$random$subject$intercept_sd <- 0.35868800362449987
  old$random$subject$slopes$cond <- 0.1195626678748333
  old$random$item$intercept_sd <- 0.2391253357496666
  old$random$item$slopes$cond <- 0.05978133393741665
  old$response$sigma <- 0.89672000906124971
  expect_equal(calibrate_response(pilotr_example("crossed_mixed_rt"), 1, tune = "all",
                                  scale_intercept = TRUE),
               old, tolerance = 1e-12)

  old <- example_spec("poisson_counts_between")
  old$fixed$intercept <- 11.377260200875366
  old$fixed$coefficients$grp <- 2.8277022992109773
  expect_equal(calibrate_response(pilotr_example("poisson_counts_between"), 2, tune = "all",
                                  scale_intercept = TRUE),
               old, tolerance = 1e-7)

  # The search still starts from a linear predictor of zero, so target 1 is still refused, in
  # words that say where the search starts.
  expect_error(
    calibrate_response(pilotr_example("poisson_counts_between"), 1, tune = "all",
                       scale_intercept = TRUE),
    paste0("with every term set to zero, the intercept included, the poisson family's residual ",
           "variance at a linear predictor of zero is 1.64493 on the latent scale, which already ",
           "meets or exceeds the target of 1; raise `target_var` above 1.64493, or hold the ",
           "intercept with scale_intercept = FALSE"),
    fixed = TRUE)

  # With an intercept of zero there is nothing to hold, and the two settings agree. A refusal
  # then has no reason to suggest holding it.
  expect_identical(calibrate_response(pilotr_example("beta_proportion"), 1, tune = "all"),
                   calibrate_response(pilotr_example("beta_proportion"), 1, tune = "all",
                                      scale_intercept = TRUE))
  expect_error(calibrate_response(pilotr_example("beta_proportion"), 0.5, tune = "all",
                                  scale_intercept = TRUE),
               "target of 0.5; raise `target_var` above 0.567646$")
})

test_that("denominator = \"random\" calibrates the total less the fixed component", {
  # The fixed effects carry 15.4% of this design's total, so the two denominators differ. The cond
  # coefficient is 0.1170 in total standard deviations, and 0.1272 in the standard deviation of
  # the random components plus the residual.
  spec <- example_spec("reading_time_continuous")
  v <- response_variance(spec)
  expect_equal(v$fixed / v$total, 0.1539, tolerance = 1e-3)
  expect_equal(calibrate_response(spec, 1, tune = "all")$fixed$coefficients$cond, 0.1170,
               tolerance = 1e-3)
  cal <- calibrate_response(spec, 1, tune = "all", denominator = "random")
  expect_equal(cal$fixed$coefficients$cond, 0.04 / sqrt(v$total - v$fixed), tolerance = 1e-10)
  expect_equal(cal$fixed$coefficients$cond, 0.1272, tolerance = 1e-3)
  expect_identical(cal$fixed$intercept, -1.2)
  vc <- response_variance(cal)
  expect_equal(vc$total - vc$fixed, 1, tolerance = 1e-10)

  # With tune = "sigma", the residual makes up what the random components leave of the target.
  sig <- calibrate_response(spec, 1, denominator = "random")
  vs <- response_variance(sig)
  expect_equal(vs$total - vs$fixed, 1, tolerance = 1e-10)
  expect_identical(sig$fixed, spec$fixed)
  expect_identical(sig$random, spec$random)
})

test_that("denominator = \"random\" makes a two-group contrast's coefficient Cohen's d", {
  skip_on_cran()   # five simulations of 40000 rows
  # 20000 subjects a group, coded -0.5 and 0.5, so the Monte Carlo standard error of a sample d
  # near 0.8 is about 0.010.
  spec <- example_spec("between_2group_gaussian")
  spec$units$subject$n <- 40000L
  spec$fixed$coefficients$grp <- 8
  se <- sqrt(2 / 20000 + 0.8^2 / (2 * 40000))

  cal <- calibrate_response(spec, 1, tune = "all", denominator = "random")
  expect_equal(cal$fixed$coefficients$grp, 0.8)
  expect_lt(abs(sample_d(cal) - cal$fixed$coefficients$grp), 3 * se)

  # The default puts the fixed component in the denominator. With tune = "sigma", a coefficient of
  # 0.8 then leaves a residual standard deviation of sqrt(1 - 0.16), and d is 0.873.
  spec$fixed$coefficients$grp <- 0.8
  total <- calibrate_response(spec, 1)
  expect_equal(total$response$sigma, sqrt(1 - 0.8^2 / 4), tolerance = 1e-4)
  expect_gt(sample_d(total) - 0.8, 3 * se)
  expect_identical(calibrate_response(spec, 1, denominator = "random")$response$sigma, 1)
})

test_that("denominator = \"random\" works for every residual kind", {
  # A residual fixed by the link: an accuracy version of the crossed design.
  acc <- example_spec("crossed_mixed_rt")
  acc$fixed$intercept <- 1
  acc$response <- list(family = "bernoulli", name = "correct")
  cal <- calibrate_response(acc, 4, tune = "all", denominator = "random")
  v <- response_variance(cal)
  expect_equal(v$total - v$fixed, 4, tolerance = 1e-10)
  expect_identical(cal$fixed$intercept, 1)
  # Without random effects, nothing in that denominator can move.
  expect_error(calibrate_response(pilotr_example("ordinal_likert_between"), 4, tune = "all",
                                  denominator = "random"),
               "this design has no random-effect variance to rescale", fixed = TRUE)
  # A residual that moves with the linear predictor, which here is all the denominator holds.
  for (nm in c("poisson_counts_between", "beta_proportion")) {
    cal <- calibrate_response(pilotr_example(nm), 1, tune = "all", denominator = "random")
    v <- response_variance(cal)
    expect_equal(v$total - v$fixed, 1, tolerance = 1e-6, info = nm)
    expect_s3_class(simulate_design(cal), "data.frame")
  }
})

test_that("calibrate_response() checks scale_intercept and denominator", {
  spec <- pilotr_example("crossed_mixed_rt")
  for (bad in list(NA, "yes", c(TRUE, FALSE), 1, NULL))
    expect_error(calibrate_response(spec, 1, tune = "all", scale_intercept = bad),
                 "`scale_intercept` must be TRUE or FALSE", fixed = TRUE)
  expect_error(calibrate_response(spec, 1, denominator = "fixed"), "should be one of")
  # tune = "sigma" scales nothing but the residual, so there is no intercept to hold.
  expect_identical(calibrate_response(spec, 1, scale_intercept = TRUE),
                   calibrate_response(spec, 1))
})
