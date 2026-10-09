# power_mixed() and precision_design() fit the model each response family needs. They fitted a
# linear mixed model to every response, so for a bernoulli or poisson outcome the estimate was a
# difference on the response scale, compared with a true value on the link scale. A design
# without random effects failed every replicate. These tests hold the fitter to the family, and
# hold the warnings for the two families that have no model of their own here yet.

# The crossed accuracy design of the audit: subjects crossed with items, an intercept of 1 and an
# effect of 0.5 on the logit scale, with by-subject and by-item intercepts alone.
bernoulli_crossed <- function(n_subject = 20, n_item = 16) build_spec(list(
  name = "accuracy", seed = 1, design_kind = "within", include_items = TRUE,
  n_subject = n_subject, n_item = n_item, factor_name = "cond", lev1 = "a", lev2 = "b",
  intercept = 1, effect = 0.5, subj_int_sd = 0.5, subj_slope_sd = 0,
  item_int_sd = 0.3, item_slope_sd = 0, family = "bernoulli", resp_name = ""))

# The same crossing for counts, with an effect of 0.3 on the log scale.
poisson_crossed <- function() build_spec(list(
  name = "counts", seed = 2, design_kind = "within", include_items = TRUE,
  n_subject = 20, n_item = 16, factor_name = "cond", lev1 = "a", lev2 = "b",
  intercept = 1, effect = 0.3, subj_int_sd = 0.3, subj_slope_sd = 0,
  item_int_sd = 0.2, item_slope_sd = 0, family = "poisson", resp_name = ""))

# The two warnings, word for word. The first reads the same in the Python twin.
ORDINAL_POWER <- paste0(
  "power_mixed() fits a linear model to the ordinal response on its own scale, while the ",
  "specification's coefficients are on the logit scale, so the mean estimate and Type M are ",
  "withheld (NA) and Type S compares signs only. For a model on the link scale, use ",
  "generate_design_analysis() or brms_bridge() in the R package.")
ORDINAL_PRECISION <- paste0(
  "precision_design() fits a linear model to the ordinal response on its own scale, while the ",
  "region of practical equivalence and the specification's coefficients are on the logit ",
  "scale, so the decision probabilities are withheld (NA). For a decision on the link scale, ",
  "use generate_design_analysis() in the R package.")

# Run `expr` and return its value with the messages of the warnings it raised, muffled.
warnings_of <- function(expr) {
  msgs <- character(0)
  value <- withCallingHandlers(expr, warning = function(w) {
    msgs <<- c(msgs, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = msgs)
}

test_that("a bernoulli response is fitted by glmer, on the logit scale of its coefficients", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  out <- power_mixed(bernoulli_crossed(), n_sims = 30)
  expect_identical(out$fitter, "lme4::glmer (binomial), Wald z")
  # Over these 30 replicates the glmer estimates have a standard deviation of about 0.19, so the
  # Monte Carlo standard error of their mean is about 0.034, and 0.1 is three of them. A linear
  # mixed model of the 0/1 response estimated the difference in probability, 0.090 here, and
  # reported a Type M of 0.21, so a design that exaggerates read as one that underestimates.
  expect_lt(abs(out$mean_estimate[["effect"]] - 0.5), 0.1)
  expect_gt(out$type_m[["effect"]], 1)
})

test_that("a poisson response is fitted by glmer, on the log scale of its coefficients", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  out <- power_mixed(poisson_crossed(), n_sims = 10)
  expect_identical(out$fitter, "lme4::glmer (poisson), Wald z")
  # The estimates have a standard deviation of about 0.05, so over 10 replicates three Monte
  # Carlo standard errors come to 0.05. The linear mixed model of the counts gave 0.90.
  expect_lt(abs(out$mean_estimate[["effect"]] - 0.3), 0.05)
})

test_that("precision_design takes a bernoulli effect's interval on the logit scale", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  pr <- precision_design(bernoulli_crossed(), rope = 0.1, n_sims = 10)
  expect_identical(pr$fitter, "lme4::glmer (binomial), Wald z")
  # A logit-scale interval is about 0.74 wide on this design. The linear model's interval, about
  # 0.12, was the width of a difference in probability, set against a region on the logit scale.
  expect_gt(pr$mean_ci_width, 0.5)
})

test_that("every shipped between-subjects example has a power path", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  # None has random effects, and lmer() refused each with "No random effects terms specified in
  # formula", so every replicate failed and power was NA in all four.
  fitters <- c(between_2group_gaussian = "stats::lm, t",
               poisson_counts_between = "stats::glm (poisson), Wald z",
               ordinal_likert_between = "stats::lm, t",
               beta_proportion = "stats::lm, t")
  for (nm in names(fitters)) {
    run <- warnings_of(power_mixed(load_spec(pilotr_example(nm)), n_sims = 5))
    out <- run$value
    expect_identical(out$n_returned, 5L, info = nm)
    expect_true(all(is.finite(out$power)), info = nm)
    expect_identical(out$fitter, fitters[[nm]], info = nm)
    # Only the two families without a model of their own warn, about their scale.
    expect_identical(length(run$warnings), if (nm %in% c("ordinal_likert_between",
                                                         "beta_proportion")) 1L else 0L,
                     info = nm)
  }
  pr <- precision_design(load_spec(pilotr_example("between_2group_gaussian")), rope = 1,
                         n_sims = 5)
  expect_identical(pr$n_returned, 5L)
  expect_true(is.finite(pr$p_meaningful))
  expect_identical(pr$fitter, "stats::lm, Wald z")
})

test_that("an ordinal response warns once and withholds the mean estimate and Type M", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  run <- warnings_of(power_mixed(load_spec(pilotr_example("ordinal_likert_between")),
                                 n_sims = 5))
  expect_identical(run$warnings, ORDINAL_POWER)
  out <- run$value
  expect_true(is.na(out$mean_estimate[["grp"]]))
  expect_true(is.na(out$type_m[["grp"]]))
  # A monotone link keeps the sign, so Type S still means what it says.
  expect_identical(out$type_s[["grp"]], 0)
  expect_true(is.finite(out$power[["grp"]]))
})

test_that("an interaction's Type S is withheld too, since its sign can change with the scale", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  # A within factor crossed with a between factor, near the top of a bounded response. The cell
  # means are 0.525, 0.802, 0.802 and 0.957, so the positive logit-scale interaction of 0.4 is a
  # difference of differences of -0.122 on the response scale. The linear model estimated it
  # below zero and significant in every replicate, and Type S read 1 for a design with no sign
  # problem at all. The main effects keep their sign, and their Type S stays.
  spec <- list(
    name = "ix", seed = 11, units = list(subject = list(n = 300)),
    factors = list(
      list(name = "A", levels = c("a1", "a2"), contrasts = list(ca = c(-0.5, 0.5)),
           vary_within = "subject"),
      list(name = "B", levels = c("b1", "b2"), contrasts = list(cb = c(-0.5, 0.5)),
           between = "subject")),
    fixed = list(intercept = 1.5, coefficients = list(ca = 1.5, cb = 1.5, `ca:cb` = 0.4)),
    random = list(subject = list(intercept_sd = 0.3)),
    response = list(family = "beta", name = "y", phi = 20))
  run <- warnings_of(power_mixed(spec, n_sims = 5))
  expect_length(run$warnings, 1L)
  out <- run$value
  expect_gt(out$n_significant[["ca_cb"]], 0L)
  expect_true(is.na(out$type_s[["ca_cb"]]))
  expect_identical(unname(out$type_s[c("ca", "cb")]), c(0, 0))
})

test_that("precision_design withholds the ROPE decisions for an ordinal response", {
  skip_if_not_installed("lme4")
  run <- warnings_of(precision_design(load_spec(pilotr_example("ordinal_likert_between")),
                                      rope = 0.1, n_sims = 5))
  expect_identical(run$warnings, ORDINAL_PRECISION)
  pr <- run$value
  withheld <- c("p_meaningful", "p_meaningful_mcse", "p_meaningful_lo", "p_meaningful_hi",
                "p_equivalent", "p_equivalent_mcse", "p_equivalent_lo", "p_equivalent_hi")
  for (col in withheld) expect_true(is.na(pr[[col]]), info = col)
  # The width stays, on the response scale, the scale of the model that was fitted.
  expect_true(is.finite(pr$mean_ci_width))
  expect_identical(pr$fitter, "stats::lm, Wald z")
})

test_that("a formula the user gives is fitted as a linear model, by lm without random terms", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  # The family chooses the fitter only for the model pilotr derives. A formula of the user's own
  # is fitted as it always was, now by lm() when it has no random terms, where lmer() refused it.
  spec <- bernoulli_crossed(8, 6)
  out <- power_mixed(spec, formula = .y ~ effect, n_sims = 2)
  expect_identical(out$fitter, "stats::lm, t")
  expect_identical(out$n_returned, 2L)
  out <- power_mixed(spec, formula = .y ~ effect + (1 | subject), n_sims = 2)
  expect_identical(out$fitter, "lmerTest::lmer, Satterthwaite t")
  # lmer() also takes a formula written as text. Read as text, it showed no random terms, so lm()
  # fitted it without its random effects.
  out <- power_mixed(spec, formula = ".y ~ effect + (1 | subject)", n_sims = 2)
  expect_identical(out$fitter, "lmerTest::lmer, Satterthwaite t")
  pr <- precision_design(spec, formula = ".y ~ effect + (1 | subject)", n_sims = 2)
  expect_identical(pr$fitter, "lme4::lmer, Wald z")
})

test_that("a glm that warns is counted as a fit with a warning, and lm fits are never singular", {
  skip_if_not_installed("lme4")
  # Complete separation: glm() reports fitted probabilities of 0 or 1.
  d <- data.frame(x = c(-2, -1, -0.5, 0.5, 1, 2), .y = c(0, 0, 0, 1, 1, 1))
  f <- pilotr:::.fit_model(.y ~ x, d, family = "bernoulli", test = TRUE, auto = TRUE)
  expect_s3_class(f$fit, "glm")
  expect_false(f$singular)
  expect_gt(length(f$messages), 0L)
  expect_false(f$converged)
  g <- pilotr:::.fit_model(.y ~ x, transform(d, .y = x + c(0.1, -0.1, 0, 0.2, -0.2, 0)),
                           family = "gaussian", test = TRUE, auto = TRUE)
  expect_s3_class(g$fit, "lm")
  expect_false(g$singular)
  expect_true(g$converged)
})

test_that("a random term is found wherever lme4 looks for one", {
  has <- pilotr:::.has_bars
  expect_true(has(.y ~ x + (1 | g)))
  expect_true(has(.y ~ x + (1 + x || g) + (1 | h)))
  expect_true(has(~ (1 | g)))
  expect_false(has(.y ~ x + z + x:z))
  expect_false(has(.y ~ 1))
  expect_true(has(".y ~ x + (1 | g)"))
  expect_false(has(".y ~ x"))
})

test_that("a power sweep carries the fitter on every row", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  curve <- power_curve_mixed(load_spec(pilotr_example("between_2group_gaussian")),
                             subject_ns = c(20, 30), n_sims = 2)
  expect_identical(curve$fitter, rep("stats::lm, t", 2))
})
