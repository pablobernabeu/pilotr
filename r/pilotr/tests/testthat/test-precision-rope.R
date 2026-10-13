# precision_design() and precision_curve() take the region of practical equivalence as a
# half-width or as an interval, and refuse anything else. The decisions compared each replicate's
# interval with `rope` and `-rope`, so a vector was recycled across the replicates, each bound
# standing in turn for the half-width, and a negative value was used as given. With an even number
# of replicates, R did not even warn.

# The design of the precision_design() example: 12 subjects crossed with 12 items and an effect of
# 0.05 on the log scale. Over 30 replicates, the range form c(-0.02, 0.02) reported a p_meaningful
# of 0.200 here, where the half-width 0.02 gave 0.067.
rope_crossed <- function() build_spec(list(
  name = "pr", seed = 1, design_kind = "within", include_items = TRUE,
  n_subject = 12, n_item = 12, factor_name = "cond", lev1 = "a", lev2 = "b",
  intercept = 6, effect = 0.05, subj_int_sd = 0.12, subj_slope_sd = 0.04, subj_corr = 0.2,
  item_int_sd = 0.08, item_slope_sd = 0.02, item_corr = -0.1,
  family = "shifted_lognormal", resp_name = "", sigma = 0.3, shift = 200))

# Two groups with a residual standard deviation of 0.02, fitted by lm(). At 300 subjects or more,
# the effect's standard error is at most 0.0023 and its interval under 0.01 wide. The interval
# lies inside a region from 0.02 to 0.06, and clear of one from 0.06 to 0.1, unless its estimate
# strays more than six standard errors from 0.04.
narrow_between <- function(n_subject = 400) build_spec(list(
  name = "narrow", seed = 5, design_kind = "between", n_subject = n_subject,
  factor_name = "group", lev1 = "a", lev2 = "b", intercept = 0, effect = 0.04,
  family = "gaussian", resp_name = "", sigma = 0.02))

ROPE_REFUSAL <- paste0("`rope` must be a non-negative half-width or an increasing pair ",
                       "c(lower, upper); got ")

# A prep function that stops the run if a replicate is ever simulated, so that a refusal is shown
# to come before any replicate.
no_replicates <- function(d) stop("a replicate was simulated before `rope` was checked")

test_that("a half-width r is the region c(-r, r), and a pair is the region it names", {
  as_rope <- pilotr:::.as_rope
  expect_identical(as_rope(0.02), c(-0.02, 0.02))
  expect_identical(as_rope(0), c(0, 0))
  expect_identical(as_rope(1L), c(-1, 1))
  expect_identical(as_rope(c(lower = -0.02, upper = 0.05)), c(-0.02, 0.05))
})

test_that("an interval is meaningful clear of the region and equivalent inside it, strictly", {
  # Six intervals against the region from -0.02 to 0.05: above it, inside it, across its upper
  # bound, below it, touching its upper bound from above and touching its lower bound from inside.
  # An interval that touches a bound is neither, as before.
  lo <- c(0.06, -0.01, 0.01, -0.09, 0.05, -0.02)
  hi <- c(0.09, 0.04, 0.07, -0.03, 0.08, 0.03)
  d <- pilotr:::.rope_decisions(lo, hi, c(-0.02, 0.05))
  expect_identical(d$meaningful, c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE))
  expect_identical(d$equivalent, c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE))
})

test_that("a rope given as c(-r, r) gives the decisions of the half-width r", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  spec <- rope_crossed()
  half <- precision_design(spec, focal = c(effect = 0.05), rope = 0.02, n_sims = 6)
  pair <- precision_design(spec, focal = c(effect = 0.05), rope = c(-0.02, 0.02), n_sims = 6)
  # Recycled, the pair reported 0.333 against the half-width's 0.167 over these six replicates.
  expect_identical(pair, half)
})

test_that("an asymmetric region decides by its own bounds, in the sweep too", {
  skip_if_not_installed("lme4")
  spec <- narrow_between()
  inside <- precision_design(spec, focal = c(effect = 0.04), rope = c(0.02, 0.06), n_sims = 6)
  expect_identical(c(inside$p_equivalent, inside$p_meaningful), c(1, 0))
  above <- precision_design(spec, focal = c(effect = 0.04), rope = c(0.06, 0.1), n_sims = 6)
  expect_identical(c(above$p_equivalent, above$p_meaningful), c(0, 1))
  curve <- precision_curve(spec, focal = c(effect = 0.04), subject_ns = c(300, 400),
                           rope = c(0.02, 0.06), n_sims = 4)
  expect_identical(curve$p_equivalent, c(1, 1))
  expect_identical(curve$p_meaningful, c(0, 0))
})

test_that("any other rope is refused before a replicate is simulated", {
  skip_if_not_installed("lme4")
  spec <- narrow_between()
  # Each value with the text the refusal shows for it.
  bad <- list("-0.02" = -0.02, "c(0.06, 0.02)" = c(0.06, 0.02), "c(1, 2, 3)" = c(1, 2, 3),
              "NA" = NA, "NA_real_" = NA_real_, "c(0.02, 0.02)" = c(0.02, 0.02), "Inf" = Inf,
              "c(-Inf, 0.02)" = c(-Inf, 0.02), "\"0.02\"" = "0.02", "TRUE" = TRUE,
              "numeric(0)" = numeric(0), "NULL" = NULL)
  for (shown in names(bad)) {
    expect_error(precision_design(spec, prep = no_replicates, rope = bad[[shown]], n_sims = 2),
                 paste0(ROPE_REFUSAL, shown), fixed = TRUE, info = shown)
  }
  expect_error(precision_curve(spec, subject_ns = c(300, 400), prep = no_replicates,
                               rope = c(0.06, 0.02), n_sims = 2),
               paste0(ROPE_REFUSAL, "c(0.06, 0.02)"), fixed = TRUE)
})
