# power_design() t-tests every row, and that test holds its level only with one row per subject
# and no shared cluster. It used to accept any Gaussian design with a two-level between factor,
# and report power for a test that did not hold its level. Each refusal below is byte-identical
# to the Python twin's, which tools/parity/validate_cross.py checks.

# A two-group Gaussian design with one row per subject, the design the t-test analyses.
one_row_per_subject <- function(n = 40, effect = 0.5) {
  build_spec(list(name = "tg", seed = 11, design_kind = "between", n_subject = n,
                  factor_name = "group", lev1 = "a", lev2 = "b", intercept = 0,
                  effect = effect, family = "gaussian", resp_name = "y", sigma = 1))
}

correlated_rows <- function(what) {
  paste0("power_design() in R and power() in Python t-test every row as an independent ",
         "observation, which is valid only when each subject contributes one row and no rows ",
         "share a cluster. This design has ", what, ", so its rows are correlated and the ",
         "t-test would overstate power and understate Type M. Use power_mixed() in the R ",
         "package, which fits the model the specification implies.")
}

test_that("a design crossed with items is refused", {
  # 30 subjects by 20 items with no true effect. Analysed row by row, 105 of 200 replicates were
  # significant, against 16 for a t-test of the subject means on the same data.
  spec <- one_row_per_subject(n = 30, effect = 0)
  spec$units$item <- list(n = 20L)
  spec$random <- list(subject = list(intercept_sd = 1), item = list(intercept_sd = 0.3))
  expect_error(power_design(spec, n_sims = 2), correlated_rows("an item unit"), fixed = TRUE)
})

test_that("subjects nested in sites are refused, and the grouping factor is named", {
  # One row per subject, but the rows of a site share its effect: 0.200 under the null.
  spec <- one_row_per_subject(n = 120, effect = 0)
  spec$random <- list(site = list(intercept_sd = 0.5, over = "subject", n = 12L))
  expect_error(power_design(spec, n_sims = 2), correlated_rows("the grouping factor 'site'"),
               fixed = TRUE)
})

test_that("a within factor is refused, and named", {
  spec <- one_row_per_subject()
  spec$factors[[2]] <- list(name = "block", levels = c("x", "y"),
                            contrasts = list(blk = c(-0.5, 0.5)), vary_within = "subject")
  expect_error(power_design(spec, n_sims = 2), correlated_rows("the within factor 'block'"),
               fixed = TRUE)
})

test_that("predictors and by-subject random effects keep one row per subject", {
  spec <- one_row_per_subject()
  spec$random <- list(subject = list(intercept_sd = 0.5))
  spec$predictors <- list(list(name = "age", varies_by = "subject", mean = 0, sd = 1),
                          list(name = "noise", varies_by = "observation", mean = 0, sd = 1))
  spec$fixed$coefficients$age <- 0.2
  r <- power_design(spec, n_sims = 20)
  expect_identical(r$true_effect, 0.5)
  expect_true(r$power >= 0 && r$power <= 1)
})

test_that("designs the backend never covered keep their messages", {
  counts <- one_row_per_subject()
  counts$response <- list(family = "poisson", name = "count")
  expect_error(power_design(counts, n_sims = 2),
               "The power backend currently handles only the gaussian two-group design.",
               fixed = TRUE)
  three <- one_row_per_subject()
  three$factors[[1]]$levels <- c("a", "b", "c")
  three$factors[[1]]$contrasts <- list(effect = c(-1, 0, 1))
  expect_error(power_design(three, n_sims = 2),
               "The power backend expects exactly one 2-level between factor.", fixed = TRUE)
})

test_that("a specification with no coefficient for the contrast has a true effect of 0", {
  # `"coefficients": {}` is a valid way to write a null design. power_design() stopped with
  # "missing value where TRUE/FALSE needed".
  spec <- one_row_per_subject()
  spec$fixed$coefficients <- list()
  r <- power_design(spec, n_sims = 20)
  expect_identical(r$true_effect, 0)
  expect_true(is.nan(r$type_s))
  expect_true(is.nan(r$type_m))
  expect_true(r$power >= 0 && r$power <= 1)
})
