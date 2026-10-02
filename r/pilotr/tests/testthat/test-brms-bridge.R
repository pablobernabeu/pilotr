# brms_bridge() returns its model rather than writing it, so that assigning the result, or
# reading it from generate_design_analysis(), costs nothing in console output. Putting the code
# on screen is the print method's job. These tests hold that split: the function silent, the
# method not, and the same four fields reaching a caller either way.
#
# They also hold the model itself to what brms accepts and to the scale of the data. brms is not a
# dependency of this package, even a suggested one, so nothing here fits a model or asks brms
# anything. tools/brms/check_bridge.R runs brms's own validation over every shipped example.

# The spec from the function's own example: crossed, with correlated by-subject and by-item
# slopes, so the bridge assembles a `|` bar and an LKJ prior as well as the plain parts.
bridge_spec <- function() {
  build_spec(list(name = "d", seed = 1, design_kind = "within",
                  include_items = TRUE, n_subject = 20, n_item = 12, factor_name = "cond",
                  lev1 = "a", lev2 = "b", intercept = 6, effect = 0.05,
                  subj_int_sd = 0.12, subj_slope_sd = 0.04, subj_corr = 0.2,
                  item_int_sd = 0.08, item_slope_sd = 0.02, item_corr = -0.1,
                  family = "shifted_lognormal", resp_name = "RT", sigma = 0.3, shift = 200))
}

# The standard deviation of every normal prior whose text contains `what`.
prior_sd <- function(priors, what) {
  hit <- grep(what, priors, value = TRUE, fixed = TRUE)
  as.numeric(sub("^prior\\(normal\\(0, ([^)]+)\\).*$", "\\1", hit))
}

# Every term that carries scale multiplied by k: the intercept, the coefficients, the
# random-effect standard deviations and the residual standard deviation.
scale_design <- function(spec, k) {
  spec[["fixed"]][["intercept"]] <- spec[["fixed"]][["intercept"]] * k
  for (nm in names(spec[["fixed"]][["coefficients"]]))
    spec[["fixed"]][["coefficients"]][[nm]] <- spec[["fixed"]][["coefficients"]][[nm]] * k
  for (g in names(spec[["random"]])) {
    spec[["random"]][[g]][["intercept_sd"]] <- spec[["random"]][[g]][["intercept_sd"]] * k
    for (nm in names(spec[["random"]][[g]][["slopes"]]))
      spec[["random"]][[g]][["slopes"]][[nm]] <- spec[["random"]][[g]][["slopes"]][[nm]] * k
  }
  spec[["response"]][["sigma"]] <- spec[["response"]][["sigma"]] * k
  spec
}

# reading_time_continuous with 10 subjects and 8 items in place of 50 and 40. Its predictors,
# interactions and random slopes are unchanged, and the widths are measured on a design one
# twenty-fifth the size, which keeps the tests that need several bridges of it quick.
small_continuous <- function() {
  spec <- load_spec(pilotr_example("reading_time_continuous"))
  spec[["units"]][["subject"]][["n"]] <- 10L
  spec[["units"]][["item"]][["n"]] <- 8L
  spec
}

# The same with a Gaussian response, which puts the coefficients on the raw scale.
gaussian_continuous <- function() {
  spec <- small_continuous()
  spec[["response"]] <- list(family = "gaussian", name = "rt", sigma = 0.25)
  spec
}

test_that("assigning a bridge writes nothing and signals nothing", {
  spec <- bridge_spec()
  # expect_silent() fails on output, on a message and on a warning alike, which is the whole
  # of the contract in one assertion.
  expect_silent(bridge <- brms_bridge(spec))
  # Then each stream separately, so a failure says which one carried the output.
  expect_length(capture.output(bridge <- brms_bridge(spec)), 0L)
  expect_length(capture.output(bridge <- brms_bridge(spec), type = "message"), 0L)
})

test_that("the result is returned visibly, so a bare call still shows the model", {
  expect_true(withVisible(brms_bridge(bridge_spec()))$visible)
})

test_that("the object carries the formula, the family, the priors and the code", {
  b <- brms_bridge(bridge_spec())
  expect_s3_class(b, "pilotr_bridge")
  expect_named(b, c("formula", "family", "priors", "code"))
  expect_identical(b$formula, "RT ~ effect + (1 + effect | subject) + (1 + effect | item)")
  expect_identical(b$family, "shifted_lognormal()")
  # On the log scale, where the shifted lognormal's coefficients live, the response's total
  # variance is 0.111926. The effect contributes 0.000626 (0.05^2 times the variance of the
  # +/-0.5 column over 480 rows, 0.25 * 480 / 479), the subjects 0.0148 (0.12^2 + 0.04^2 * 0.25),
  # the items 0.0065 (0.08^2 + 0.02^2 * 0.25) and the residual 0.09. The SD is then 0.3346, and
  # the effect's prior is half of it per SD of its column, 0.5 * 0.3346 / 0.5005 = 0.3342.
  expect_identical(b$priors, c(
    'prior(normal(0, 0.3342), class = "b", coef = "effect")',
    'prior(normal(0, 0.3346), class = "sd")',
    'prior(lkj(2), class = "cor")'))
  # The code is the other three assembled, not a fifth thing that could drift from them.
  expect_true(grepl(b$formula, b$code, fixed = TRUE))
  expect_true(grepl(paste0("family = ", b$family), b$code, fixed = TRUE))
  for (p in b$priors) expect_true(grepl(p, b$code, fixed = TRUE))
})

# brms refuses a prior on a parameter the model does not contain, and a design without random
# effects has no standard deviations to put one on. The prior used to be emitted regardless, so
# brm() stopped on four of the eight shipped examples before it sampled anything. The intercept
# prior is left to brms, whose default is centred on the data. A unit-scale one on a response
# near 100 pulled the intercept to about 4, and the residual SD up to about 95 to make up for it.
test_that("every example has an SD prior only with random effects, and no intercept prior", {
  for (ex in pilotr_example()) {
    spec <- load_spec(pilotr_example(ex))
    b <- brms_bridge(spec)
    expect_identical(sum(grepl('class = "sd"', b$priors, fixed = TRUE)),
                     if (length(spec[["random"]])) 1L else 0L, info = ex)
    expect_false(any(grepl('class = "Intercept"', b$priors, fixed = TRUE)), info = ex)
  }
  b <- brms_bridge(load_spec(pilotr_example("between_2group_gaussian")))
  expect_false(grepl('class = "sd"', b$code, fixed = TRUE))
})

test_that("the emitted call says that its data come from model_data()", {
  b <- brms_bridge(load_spec(pilotr_example("between_2group_gaussian")))
  expect_true(grepl("data   = your_data,   # model_data(spec, simulate_design(spec))", b$code,
                    fixed = TRUE))
  # The comment sits inside the call, so the call still parses.
  expect_type(parse(text = b$code), "expression")
})

# The grp column is -0.5 for 32 subjects and 0.5 for the other 32, so its variance is
# 0.25 * 64 / 63 and the effect of 5 adds 25 times that to the residual variance of 100. Half of
# the response's SD per SD of the column is then sqrt(25 * 0.25 + 100 * 63 / 64) = 10.23 points.
# The old unit-scale N(0, 0.5) held the effect near zero whatever the data said, which made the
# Savage-Dickey Bayes factor about 1.
test_that("a Gaussian effect's prior is on the response's own scale", {
  b <- brms_bridge(load_spec(pilotr_example("between_2group_gaussian")))
  expect_identical(b$priors, 'prior(normal(0, 10.23), class = "b", coef = "grp")')
})

test_that("multiplying a Gaussian design's scale by k multiplies every prior width by k", {
  designs <- list(nested_clusters = load_spec(pilotr_example("nested_clusters")),
                  partial_crossing = load_spec(pilotr_example("partial_crossing")),
                  continuous = gaussian_continuous())
  for (nm in names(designs)) {
    spec <- designs[[nm]]
    one <- brms_bridge(spec)$priors
    ten <- brms_bridge(scale_design(spec, 10))$priors
    for (what in c('class = "b"', 'class = "sd"')) {
      w1 <- prior_sd(one, what)
      expect_gt(length(w1), 0L)
      expect_equal(prior_sd(ten, what), 10 * w1, tolerance = 1e-3, info = paste(nm, what))
    }
  }
})

# brms names an interaction after its parts in the order they first appear in the formula, as
# stats::terms() does. reading_time_continuous lists age before cond, so its "cond:age" is the
# coefficient brms calls "age:cond", and a prior written for "cond:age" stopped brm().
test_that("an interaction is written as brms names it, in the formula and in its prior", {
  b <- brms_bridge(load_spec(pilotr_example("reading_time_continuous")))
  expect_true(grepl(" + age:cond + ", b$formula, fixed = TRUE))
  expect_false(grepl("cond:age", b$formula, fixed = TRUE))
  expect_length(grep('coef = "age:cond"', b$priors, fixed = TRUE), 1L)
  expect_length(grep("cond:age", b$priors, fixed = TRUE), 0L)
  # Interactions whose parts already appear in that order keep their names.
  expect_length(grep('coef = "SyntaxPC:age"', b$priors, fixed = TRUE), 1L)
  expect_length(grep('coef = "CoherencePC:age"', b$priors, fixed = TRUE), 1L)
})

# A specification may give one interaction two keys, and the simulation adds their effects on the
# one product column. brms has a single coefficient for them, and it refused the code the bridge
# wrote, with two priors on that coefficient ("Duplicated prior specifications are not allowed").
test_that("one interaction under two keys is refused, since brms gives it one coefficient", {
  spec <- small_continuous()
  spec[["fixed"]][["coefficients"]][["age:cond"]] <- 0.01
  msg <- paste("the keys 'cond:age' and 'age:cond' name one interaction, which brms estimates",
               "as the single coefficient 'age:cond'; write each interaction once")
  expect_error(brms_bridge(spec), msg, fixed = TRUE)
  # The design analysis takes its model from the bridge, so it stops on the same grounds.
  expect_error(generate_design_analysis(spec, focal = "cond"), msg, fixed = TRUE)
})

test_that("prior_scale and interaction_scale set the standardised widths", {
  spec <- small_continuous()   # has interaction terms
  main <- c('coef = "cond"', 'coef = "age"')
  inter <- c('coef = "age:cond"', 'coef = "SyntaxPC:age"')
  base <- brms_bridge(spec)$priors
  wide <- brms_bridge(spec, prior_scale = 1)$priors
  tight <- brms_bridge(spec, prior_scale = 1, interaction_scale = 0.1)$priors
  for (w in main) {
    expect_equal(prior_sd(wide, w), 2 * prior_sd(base, w), tolerance = 1e-3, info = w)
    expect_equal(prior_sd(tight, w), 2 * prior_sd(base, w), tolerance = 1e-3, info = w)
  }
  for (w in inter) {
    # An unset interaction_scale is still half the main-effect scale.
    expect_equal(prior_sd(wide, w), 2 * prior_sd(base, w), tolerance = 1e-3, info = w)
    expect_equal(prior_sd(tight, w), 0.4 * prior_sd(base, w), tolerance = 1e-3, info = w)
  }
})

# For a link family, the coefficients are on the link scale, which does not depend on the
# response's units, so the widths are prior_scale and interaction_scale as given.
test_that("a link family keeps unit-scale priors on its coefficients and random SDs", {
  spec <- small_continuous()
  spec[["response"]] <- list(family = "bernoulli", name = "correct")
  wide <- brms_bridge(spec, prior_scale = 1)
  expect_true('prior(normal(0, 1), class = "b", coef = "cond")' %in% wide$priors)
  expect_true('prior(normal(0, 0.5), class = "b", coef = "age:cond")' %in% wide$priors)
  expect_true('prior(normal(0, 1), class = "sd")' %in% wide$priors)
  tight <- brms_bridge(spec, prior_scale = 1, interaction_scale = 0.1)
  expect_true('prior(normal(0, 0.1), class = "b", coef = "age:cond")' %in% tight$priors)
  expect_true('prior(normal(0, 1), class = "b", coef = "cond")' %in% tight$priors)
  expect_identical(brms_bridge(load_spec(pilotr_example("poisson_counts_between")))$priors,
                   'prior(normal(0, 0.5), class = "b", coef = "grp")')
})

# A column that does not vary has no standard deviation to divide by. Its coefficient cannot be
# told apart from the intercept, and an infinite prior width would not be valid code.
test_that("a column that does not vary is reported and given a finite prior", {
  spec <- load_spec(pilotr_example("between_2group_gaussian"))
  spec[["predictors"]] <- list(list(name = "x", varies_by = "subject", mean = 2, sd = 0))
  spec[["fixed"]][["coefficients"]][["x"]] <- 0.3
  expect_warning(b <- brms_bridge(spec), "'x' does not vary in this design")
  expect_equal(prior_sd(b$priors, 'coef = "x"'),
               0.5 * sqrt(response_variance(spec)$total), tolerance = 1e-3)
})

# The scales now multiply measured standard deviations, so a scale that is not a positive number
# would be written into the code as a width brms cannot use.
test_that("the scales and a specification given as a list are checked before anything is emitted", {
  spec <- load_spec(pilotr_example("poisson_counts_between"))
  expect_error(brms_bridge(spec, prior_scale = 0), "`prior_scale` must be a single positive")
  expect_error(brms_bridge(spec, prior_scale = c(0.5, 1)), "`prior_scale` must be")
  expect_error(brms_bridge(spec, interaction_scale = -1), "`interaction_scale` must be")
  # An integer scale is a number like any other.
  expect_identical(brms_bridge(spec, prior_scale = 1L)$priors,
                   'prior(normal(0, 1), class = "b", coef = "grp")')
  # Its priors are measured on the data a specification simulates, so a list is validated, as a
  # path always was.
  bad <- spec
  bad[["fixed"]][["coefficients"]] <- list(gpr = 0.4)
  expect_error(brms_bridge(bad), "invalid design specification")
})

# A design with no coefficients and no random effects is an intercept-only model, which needs a
# right-hand side of 1 and leaves brms no prior to be told about.
test_that("an intercept-only design emits `~ 1` and no prior argument", {
  spec <- load_spec(pilotr_example("between_2group_gaussian"))
  spec[["fixed"]][["coefficients"]] <- structure(list(), names = character(0))
  b <- brms_bridge(spec)
  expect_identical(b$formula, "score ~ 1")
  expect_length(b$priors, 0L)
  expect_false(grepl("prior", b$code, fixed = TRUE))
  expect_type(parse(text = b$code), "expression")
})

test_that("printing writes the code to standard output as one block", {
  b <- brms_bridge(bridge_spec())
  out <- capture.output(print(b))
  # cat(code, "\n") writes the code, a separating space and a newline, so the capture is the
  # code exactly, bar that trailing space on the last line.
  expect_identical(paste(out, collapse = "\n"), paste0(b$code, " "))
  # Nothing on the message stream, which knitr collects separately: a method split across the
  # two streams renders one printed object as two boxes on the documentation site. Standard
  # output is sunk here so that the code itself does not clutter the reporter.
  con <- file(nullfile(), open = "wt")
  sink(con)
  on.exit({ sink(); close(con) }, add = TRUE)
  expect_length(capture.output(print(b), type = "message"), 0L)
})

test_that("printing returns its input invisibly", {
  b <- brms_bridge(bridge_spec())
  con <- file(nullfile(), open = "wt")
  sink(con)
  on.exit({ sink(); close(con) }, add = TRUE)
  expect_invisible(print(b))
  expect_identical(print(b), b)
})
