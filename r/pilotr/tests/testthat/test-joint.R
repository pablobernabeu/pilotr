# A set of coefficients is tested jointly. A factor with three levels has two contrast columns,
# and its interaction with a two-level factor has two product columns, so the effect a study
# preregisters is often a test of two coefficients at once, whose power differs from that of
# either coefficient alone. `focal` took coefficient names one by one and tested each by itself.
# A named list of sets now plans the joint test: Satterthwaite's F from lmerTest::contest() for a
# linear mixed model, a Wald chi-square from the coefficients' covariance for glmer() and glm(),
# and the F test for lm().

# A 3 x 2 within design, crossed with items: A has three levels, coded by the two contrast
# columns A1 and A2, and B two. The interaction is the pair A1_B and A2_B.
three_by_two <- function(seed = 7, n_subject = 12, n_item = 6) list(
  name = "three_by_two", seed = seed,
  units = list(subject = list(n = n_subject), item = list(n = n_item)),
  factors = list(
    list(name = "A", levels = c("a1", "a2", "a3"),
         contrasts = list(A1 = c(-1/3, 2/3, -1/3), A2 = c(-1/3, -1/3, 2/3)),
         vary_within = c("subject", "item")),
    list(name = "B", levels = c("b1", "b2"), contrasts = list(B = c(-0.5, 0.5)),
         vary_within = c("subject", "item"))),
  fixed = list(intercept = 0, coefficients = list(A1 = 0.2, A2 = 0.1, B = 0.2,
                                                  "A1:B" = 0.4, "A2:B" = 0.1)),
  random = list(subject = list(intercept_sd = 0.4), item = list(intercept_sd = 0.3)),
  response = list(family = "gaussian", name = "y", sigma = 1))

# Three groups between subjects, with no random terms, so lm() fits it, or glm() for counts.
three_groups <- function(family = "gaussian") list(
  name = "three_groups", seed = 3, units = list(subject = list(n = 36)),
  factors = list(list(name = "group", levels = c("g1", "g2", "g3"),
                      contrasts = list(G2 = c(0, 1, 0), G3 = c(0, 0, 1)),
                      between = "subject")),
  fixed = list(intercept = if (family == "poisson") 1 else 0,
               coefficients = list(G2 = 0.3, G3 = 0.5)),
  random = list(),
  response = if (family == "poisson") list(family = "poisson", name = "y")
             else list(family = "gaussian", name = "y", sigma = 1))

# Three conditions within subjects and items, for an accuracy response fitted by glmer().
three_bernoulli <- function() {
  spec <- three_by_two(seed = 2, n_subject = 12, n_item = 8)
  spec$factors <- spec$factors[1]
  spec$fixed <- list(intercept = 1, coefficients = list(A1 = 0.5, A2 = -0.3))
  spec$response <- list(family = "bernoulli", name = "y")
  spec
}

# The first replicate's data and fit, as power_mixed() makes them.
first_fit <- function(spec) {
  s <- spec
  s$seed <- replicate_seeds(spec$seed, 1)[1]
  d <- model_data(spec, simulate_design(s))
  f <- pilotr:::.fit_model(model_formula(spec), d, family = spec$response$family, test = TRUE)
  list(data = d, fit = f$fit)
}

# One replicate of the loop behind power_mixed(), with the sets given or with each coefficient
# tested by itself.
replicate_record <- function(spec, fnames, sets = NULL) {
  seeds <- replicate_seeds(spec$seed, 1)
  pilotr:::.design_rep(1L, spec, seeds, prep = pilotr:::.default_prep(spec),
                       formula = model_formula(spec), fnames = fnames, test = TRUE,
                       auto = TRUE, sets = sets)
}

test_that("a set of one coefficient is the single-coefficient test in every fitter", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  designs <- list(
    "lmerTest::lmer" = list(spec = three_by_two(), coef = "B", stat = "t value"),
    "stats::lm" = list(spec = three_groups(), coef = "G3", stat = "t value"),
    "stats::glm" = list(spec = three_groups("poisson"), coef = "G3", stat = "z value"),
    "lme4::glmer" = list(spec = three_bernoulli(), coef = "A1", stat = "z value"))
  for (fitter in names(designs)) {
    x <- designs[[fitter]]
    fit <- first_fit(x$spec)$fit
    co <- summary(fit)$coefficients
    pcol <- intersect(c("Pr(>|t|)", "Pr(>|z|)"), colnames(co))
    joint <- pilotr:::.joint_test(fit, x$coef)
    # F is t squared, a chi-square on one degree of freedom is z squared, and the p-values agree.
    expect_equal(joint$statistic, co[x$coef, x$stat]^2, info = fitter)
    expect_equal(joint$df1, 1, info = fitter)
    expect_equal(joint$p, co[x$coef, pcol], info = fitter)
    # Through the replicate loop, the set's p-value is the coefficient's.
    one <- replicate_record(x$spec, x$coef)
    set <- replicate_record(x$spec, "only", sets = list(only = x$coef))
    expect_equal(set$p[["only"]], one$p[[x$coef]], info = fitter)
  }
})

test_that("the joint test of a linear mixed model is lmerTest's Satterthwaite F", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  spec <- three_by_two()
  fit <- first_fit(spec)$fit
  fe <- lme4::fixef(fit)
  L <- matrix(0, 2, length(fe))
  L[cbind(1:2, match(c("A1_B", "A2_B"), names(fe)))] <- 1
  ref <- lmerTest::contest(fit, L, joint = TRUE)
  joint <- pilotr:::.joint_test(fit, c("A1_B", "A2_B"))
  expect_identical(joint$statistic, ref[["F value"]])
  expect_identical(c(joint$df1, joint$df2), c(ref[["NumDF"]], ref[["DenDF"]]))
  expect_identical(joint$p, ref[["Pr(>F)"]])
  rec <- replicate_record(spec, c("interaction", "A"),
                          sets = list(interaction = c("A1_B", "A2_B"), A = c("A1", "A2")))
  expect_identical(rec$p[["interaction"]], ref[["Pr(>F)"]])
  # A set has no single estimate.
  expect_true(all(is.na(c(rec$est, rec$se))))
  expect_identical(rec$present, c(interaction = TRUE, A = TRUE))
})

test_that("the joint test of lm() is the F test of the nested models", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  d <- first_fit(three_groups())$data
  full <- stats::lm(.y ~ G2 + G3, data = d)
  ref <- stats::anova(stats::lm(.y ~ 1, data = d), full)
  joint <- pilotr:::.joint_test(full, c("G2", "G3"))
  expect_equal(joint$statistic, ref[["F"]][2])
  expect_identical(c(joint$df1, joint$df2), c(2, stats::df.residual(full)))
  expect_equal(joint$p, ref[["Pr(>F)"]][2])
})

test_that("the joint test of glm() and glmer() is the Wald chi-square of the coefficients", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  for (spec in list(three_groups("poisson"), three_bernoulli())) {
    fit <- first_fit(spec)$fit
    set <- setdiff(pilotr:::.default_focal(spec), "(Intercept)")
    b <- if (inherits(fit, "merMod")) lme4::fixef(fit)[set] else stats::coef(fit)[set]
    V <- as.matrix(stats::vcov(fit))[set, set]
    wald <- drop(t(b) %*% solve(V) %*% b)
    joint <- pilotr:::.joint_test(fit, set)
    expect_equal(joint$statistic, wald, info = class(fit)[1])
    expect_identical(joint$df1, 2, info = class(fit)[1])
    expect_true(is.na(joint$df2), info = class(fit)[1])
    expect_equal(joint$p, stats::pchisq(wald, 2, lower.tail = FALSE), info = class(fit)[1])
  }
})

test_that("power_mixed() plans the joint tests of a 3 x 2 design", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  spec <- three_by_two()
  sets <- list(interaction = c("A1_B", "A2_B"), A = c("A1", "A2"))
  out <- power_mixed(spec, focal = sets, n_sims = 6)
  expect_identical(out$fitter, "lmerTest::lmer, Satterthwaite F")
  expect_identical(out$sets, sets)
  expect_identical(names(out$power), c("interaction", "A"))
  expect_identical(out$n_returned, 6L)
  expect_true(all(out$power >= 0 & out$power <= 1))
  expect_identical(out$power, out$n_significant / 6)
  # Each power carries its Monte Carlo error and Wilson interval.
  rate <- pilotr:::.rate_with_error(out$n_significant[["interaction"]], 6L, "power")
  expect_identical(out$power_mcse[["interaction"]], rate$power_mcse)
  expect_identical(c(out$power_lo[["interaction"]], out$power_hi[["interaction"]]),
                   c(rate$power_lo, rate$power_hi))
  # A joint test has no single estimate or direction, so the design-analysis quantities are NA.
  expect_true(all(is.na(c(out$true_effect, out$mean_estimate, out$type_s, out$type_s_lo,
                          out$type_s_hi, out$type_m, out$type_m_mcse))))
  # The fits are those of the default analysis.
  single <- power_mixed(spec, focal = c("A1_B", "A2_B"), n_sims = 6)
  expect_identical(out[c("n_returned", "n_converged", "n_singular", "n_warning")],
                   single[c("n_returned", "n_converged", "n_singular", "n_warning")])
  expect_null(single$sets)
})

test_that("a set of one coefficient has the power of that coefficient", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  spec <- three_by_two()
  one <- power_mixed(spec, focal = "A1_B", n_sims = 6)
  set <- power_mixed(spec, focal = list(A1_B = "A1_B"), n_sims = 6)
  expect_identical(set$n_significant, one$n_significant)
  expect_identical(set$power, one$power)
  expect_identical(set$fitter, "lmerTest::lmer, Satterthwaite F")
})

test_that("sets are tested on the fit a simplification rule keeps", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  spec <- three_by_two()
  ladder <- c("maximal", "intercepts")
  sets <- list(interaction = c("A1_B", "A2_B"), A1_B = "A1_B")
  out <- power_mixed(spec, focal = sets, formula = ladder, n_sims = 4)
  expect_identical(out$sets, sets)
  expect_identical(out$fitter, "lmerTest::lmer, Satterthwaite F")
  expect_equal(sum(out$rung_shares), 1)
  # A set of one coefficient keeps the power of that coefficient under the same rule.
  one <- power_mixed(spec, focal = "A1_B", formula = ladder, n_sims = 4)
  expect_identical(out$n_significant[["A1_B"]], one$n_significant[["A1_B"]])
  expect_identical(out$rung_shares, one$rung_shares)
  # The sweep frame names each set's coefficients and gives each rung's share.
  frame <- pilotr:::.as_power_frame(out)
  expect_identical(names(frame)[1:2], c("effect", "coefficients"))
  expect_identical(tail(names(frame), 2), c("rung_maximal", "rung_intercepts"))
})

test_that("a sweep over sets has one row per set and point, naming the coefficients", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  sets <- list(group = c("G2", "G3"), G3 = "G3")
  curve <- power_curve_mixed(three_groups(), subject_ns = c(30, 36), focal = sets, n_sims = 4)
  expect_identical(curve$effect, rep(c("group", "G3"), 2))
  expect_identical(curve$coefficients, rep(c("G2, G3", "G3"), 2))
  expect_identical(curve$fitter, rep("stats::lm, F", 4))
  expect_true(all(is.na(curve$type_s)))
  # Without sets, the frame keeps its columns.
  plain <- power_curve_mixed(three_groups(), subject_ns = 30, n_sims = 2)
  expect_false("coefficients" %in% names(plain))
  expect_identical(plain$fitter, rep("stats::lm, t", 2))
})

test_that("a set naming a coefficient the model lacks is reported", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  expect_warning(
    out <- power_mixed(three_groups(), focal = list(group = c("G2", "G4")), n_sims = 2),
    "focal set 'group' was never tested, since the fitted model lacks 'G4'", fixed = TRUE)
  expect_true(is.na(out$power[["group"]]))
})

test_that("malformed sets are refused before any fitting", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  spec <- three_groups()
  for (bad in list(list(c("G2", "G3")), list(a = c("G2", "G3"), c("G2")),
                   list(a = "G2", a = "G3"), structure(list(), names = character(0))))
    expect_error(power_mixed(spec, focal = bad, n_sims = 2),
                 "every set needs its own name", fixed = TRUE)
  for (bad in list(list(a = character(0)), list(a = 1:2), list(a = c("G2", NA)),
                   list(a = c("G2", "G2")), list(a = c("G2", ""))))
    expect_error(power_mixed(spec, focal = bad, n_sims = 2),
                 "set 'a' of `focal` must name one or more distinct coefficients", fixed = TRUE)
})

test_that("a joint test is two-sided, and precision_design() refuses sets", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  spec <- three_groups()
  sets <- list(group = c("G2", "G3"))
  msg <- "a set of coefficients is tested jointly"
  expect_error(power_mixed(spec, focal = sets, n_sims = 2, alternative = "greater"),
               msg, fixed = TRUE)
  expect_error(power_curve_mixed(spec, subject_ns = 30, focal = sets, n_sims = 2,
                                 alternative = "less"), msg, fixed = TRUE)
  expect_error(precision_design(spec, focal = sets, rope = 0.1, n_sims = 2),
               "precision_design() builds an interval around each coefficient", fixed = TRUE)
  expect_error(precision_curve(spec, focal = sets, subject_ns = 30, rope = 0.1, n_sims = 2),
               "precision_curve() builds an interval around each coefficient", fixed = TRUE)
})

# A result with sets, built by hand as in test-print.R, since the print method only reads it.
fake_joint <- function() {
  out <- list(
    n_sims = 200L, alpha = 0.05, fitter = "lmerTest::lmer, Satterthwaite F",
    n_attempted = 200L, n_returned = 200L, n_converged = 196L, n_singular = 4L, n_warning = 0L,
    power = c(interaction = 0.46, A = 0.66), power_mcse = c(interaction = 0.035, A = 0.033),
    power_lo = c(interaction = 0.392, A = 0.592), power_hi = c(interaction = 0.529, A = 0.722),
    n_significant = c(interaction = 92L, A = 132L),
    true_effect = c(interaction = NA_real_, A = NA_real_),
    mean_estimate = c(interaction = NA_real_, A = NA_real_),
    type_s = c(interaction = NA_real_, A = NA_real_),
    type_s_lo = c(interaction = NA_real_, A = NA_real_),
    type_s_hi = c(interaction = NA_real_, A = NA_real_),
    type_m = c(interaction = NA_real_, A = NA_real_),
    type_m_mcse = c(interaction = NA_real_, A = NA_real_),
    sets = list(interaction = c("A1_B", "A2_B"), A = c("A1", "A2")))
  class(out) <- c("pilotr_power", "list")
  out
}

test_that("a result with sets prints each set's coefficients and no Type S or Type M table", {
  out <- capture.output(print(fake_joint()))
  heading <- grep("each set of coefficients tested jointly", out, fixed = TRUE)
  expect_length(heading, 1L)
  expect_true(any(grepl("^\\s+interaction: A1_B, A2_B$", out)))
  expect_true(any(grepl("^\\s+A: A1, A2$", out)))
  expect_true(any(grepl("[0.392, 0.529]", out, fixed = TRUE)))
  expect_false(any(grepl("type_s", out, fixed = TRUE)))
  expect_length(grep("Type S and Type M are not defined for a joint test", out, fixed = TRUE), 1L)
  frame <- pilotr:::.as_power_frame(fake_joint())
  expect_identical(frame$effect, c("interaction", "A"))
  expect_identical(frame$coefficients, c("A1_B, A2_B", "A1, A2"))
  expect_identical(names(frame)[1:2], c("effect", "coefficients"))
})
