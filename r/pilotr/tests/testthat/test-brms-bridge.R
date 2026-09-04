# brms_bridge() returns its model rather than writing it, so that assigning the result, or
# reading it from generate_design_analysis(), costs nothing in console output. Putting the code
# on screen is the print method's job. These tests hold that split: the function silent, the
# method not, and the same four fields reaching a caller either way.

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
  expect_identical(b$priors, c(
    'prior(normal(0, 2.5), class = "Intercept")',
    'prior(normal(0, 0.5), class = "b", coef = "effect")',
    'prior(normal(0, 1), class = "sd")',
    'prior(lkj(2), class = "cor")'))
  # The code is the other three assembled, not a fifth thing that could drift from them.
  expect_true(grepl(b$formula, b$code, fixed = TRUE))
  expect_true(grepl(paste0("family = ", b$family), b$code, fixed = TRUE))
  for (p in b$priors) expect_true(grepl(p, b$code, fixed = TRUE))
})

test_that("prior_scale and interaction_scale still set the prior widths", {
  spec <- load_spec(pilotr_example("reading_time_continuous"))   # has interaction terms
  wide <- brms_bridge(spec, prior_scale = 1)
  expect_true('prior(normal(0, 1), class = "b", coef = "cond")' %in% wide$priors)
  # An unset interaction_scale is still half the main-effect scale.
  expect_true('prior(normal(0, 0.5), class = "b", coef = "cond:age")' %in% wide$priors)
  tight <- brms_bridge(spec, prior_scale = 1, interaction_scale = 0.1)
  expect_true('prior(normal(0, 0.1), class = "b", coef = "cond:age")' %in% tight$priors)
  expect_true('prior(normal(0, 1), class = "b", coef = "cond")' %in% tight$priors)
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
