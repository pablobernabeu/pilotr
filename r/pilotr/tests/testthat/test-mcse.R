# The Monte Carlo error of every rate the power functions report, held bit for bit against the
# Python twin.
#
# mcse.R and python/pilotr/mcse.py compute each quantity with the same operations in the same
# order, so for the same counts they return the same doubles. The pins below are those doubles,
# written as hexadecimal floating-point strings in the form Python's float.hex() gives them, and
# python/tests/test_mcse.py holds the same strings and checks them against this file. They are
# compared here by value: as.numeric() reads a hexadecimal string exactly, whereas the text that
# sprintf("%a") writes differs between platforms and from Python's (0x1p+0 against
# 0x1.0000000000000p+0 for 1).

rate_pins <- utils::read.table(
  header = TRUE, colClasses = c("integer", "integer", rep("character", 4)), text = "
count  n  p                     mcse                  lo                    hi
    0 60  0x0.0p+0              0x0.0p+0              0x0.0p+0              0x1.eced8522e178fp-5
    1 60  0x1.1111111111111p-6  0x1.0ec813a58caffp-6  0x1.826bd3cdd1810p-9  0x1.6ab4c3dd18080p-4
    8 50  0x1.47ae147ae147bp-3  0x1.a8b8c95a80272p-5  0x1.55803137393f1p-4  0x1.23fc4e7a78f78p-2
   16 20  0x1.999999999999ap-1  0x1.6e5b7d16657e0p-4  0x1.2affc35bd30ebp-1  0x1.d6b409f2188edp-1
   30 60  0x1.0000000000000p-1  0x1.08654a2d4f6dap-4  0x1.82681a193ffe0p-2  0x1.3ecbf2f360010p-1
   59 60  0x1.f777777777777p-1  0x1.0ec813a58cb06p-6  0x1.d2a967845cfefp-1  0x1.fe7d942c322e7p-1
   60 60  0x1.0000000000000p+0  0x0.0p+0              0x1.e13127add1e86p-1  0x1.0000000000000p+0
")

# Seven exaggeration ratios, and the standard error and Type S interval the twins pin for them.
type_m_ratios <- c(1.82, 2.41, 1.37, 3.05, 1.96, 2.73, 1.51)
type_m_pin <- "0x1.e6a0d0061eecfp-3"
type_s_pin <- c("0x1.6f49c628ba020p-6", "0x1.e23082fab5585p-2")   # 1 wrong sign in 8

test_that("each rate, its standard error and its Wilson bounds match the pinned doubles", {
  for (i in seq_len(nrow(rate_pins))) {
    pin <- rate_pins[i, ]
    got <- pilotr:::.rate_with_error(pin$count, pin$n, "p")
    expect_identical(names(got), c("p", "p_mcse", "p_lo", "p_hi"))
    expect_identical(unname(unlist(got)), as.numeric(unlist(pin[c("p", "mcse", "lo", "hi")])),
                     info = sprintf("%d of %d", pin$count, pin$n))
  }
})

test_that("type_m_mcse is the standard deviation of the ratios over the root of their number", {
  expect_identical(pilotr:::.type_m_mcse(type_m_ratios), as.numeric(type_m_pin))
  # sd() accumulates in long double where the platform has one, so it can differ from the plain
  # double fold in the last bits, and it agrees to within that.
  expect_equal(pilotr:::.type_m_mcse(type_m_ratios),
               stats::sd(type_m_ratios) / sqrt(length(type_m_ratios)), tolerance = 1e-14)
  # One ratio has no standard deviation.
  expect_identical(pilotr:::.type_m_mcse(1.4), NA_real_)
  expect_identical(pilotr:::.type_m_mcse(numeric(0)), NA_real_)
})

test_that("the Type S interval is a Wilson interval, withheld without a true effect", {
  expect_identical(pilotr:::.type_s_interval(1L, 8L, 0.4), as.numeric(type_s_pin))
  expect_identical(pilotr:::.type_s_interval(1L, 8L, -0.4), as.numeric(type_s_pin))
  expect_identical(pilotr:::.type_s_interval(0L, 8L, 0), c(NA_real_, NA_real_))
  expect_identical(pilotr:::.type_s_interval(0L, 8L, NA_real_), c(NA_real_, NA_real_))
  expect_identical(pilotr:::.type_s_interval(0L, 0L, 0.4), c(NA_real_, NA_real_))
})

two_group <- function(effect = 5) {
  build_spec(list(name = "mc", seed = 2024, design_kind = "between", n_subject = 64,
                  factor_name = "group", lev1 = "control", lev2 = "treatment", intercept = 100,
                  effect = effect, family = "gaussian", resp_name = "", sigma = 10))
}

test_that("power_design() reports the Monte Carlo error of power, Type S and Type M", {
  r <- power_design(two_group(), n_sims = 100)
  expect_identical(names(r), c("n_sims", "alpha", "power", "power_mcse", "power_lo", "power_hi",
                               "n_significant", "true_effect", "mean_estimate", "type_s",
                               "type_s_lo", "type_s_hi", "type_m", "type_m_mcse"))
  # Still a plain list, with no class or print method, so code that reads it is unaffected.
  expect_identical(class(r), "list")
  expect_identical(r[c("power", "power_mcse", "power_lo", "power_hi")],
                   pilotr:::.rate_with_error(r$n_significant, r$n_sims, "power"))
  wrong <- round(r$type_s * r$n_significant)
  expect_identical(c(r$type_s_lo, r$type_s_hi),
                   pilotr:::.type_s_interval(wrong, r$n_significant, r$true_effect))
  expect_true(r$type_m_mcse > 0 && r$type_m_mcse < r$type_m)
})

test_that("Type M's standard error needs two significant replicates, and Type S's interval an effect", {
  # With alpha = 1, every replicate is significant, so n_sims sets the count exactly.
  one <- power_design(two_group(), n_sims = 1, alpha = 1)
  expect_identical(one$n_significant, 1L)
  expect_false(is.nan(one$type_m))
  expect_true(is.nan(one$type_m_mcse))
  expect_false(anyNA(c(one$type_s_lo, one$type_s_hi)))
  two <- power_design(two_group(), n_sims = 2, alpha = 1)
  expect_identical(two$n_significant, 2L)
  expect_true(is.finite(two$type_m_mcse))
  # Type S and Type M are undefined at a zero effect, and so is their uncertainty. Power keeps its
  # own, since it is defined there.
  null <- power_design(two_group(effect = 0), n_sims = 50)
  expect_true(all(is.nan(c(null$type_s, null$type_s_lo, null$type_s_hi, null$type_m,
                           null$type_m_mcse))))
  expect_false(anyNA(c(null$power_mcse, null$power_lo, null$power_hi)))
})

test_that("power_mixed() gives every effect a Type S interval and a Type M standard error", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  spec <- build_spec(list(name = "mc", seed = 1, design_kind = "within",
                          include_items = TRUE, n_subject = 10, n_item = 8,
                          factor_name = "cond", lev1 = "a", lev2 = "b",
                          intercept = 6, effect = 0.3, subj_int_sd = 0.12,
                          subj_slope_sd = 0, item_int_sd = 0.08, item_slope_sd = 0,
                          family = "gaussian", resp_name = "", sigma = 0.3))
  out <- power_mixed(spec, n_sims = 6)
  expect_true(out$n_significant[["effect"]] >= 2L)
  wrong <- round(out$type_s[["effect"]] * out$n_significant[["effect"]])
  expect_identical(unname(c(out$type_s_lo, out$type_s_hi)),
                   pilotr:::.type_s_interval(wrong, out$n_significant[["effect"]],
                                             out$true_effect[["effect"]]))
  expect_identical(names(out$type_m_mcse), "effect")
  expect_true(out$type_m_mcse[["effect"]] > 0)
  # The sweeps carry the same quantities as columns.
  frame <- pilotr:::.as_power_frame(out)
  expect_identical(frame$type_s_lo, unname(out$type_s_lo))
  expect_identical(frame$type_s_hi, unname(out$type_s_hi))
  expect_identical(frame$type_m_mcse, unname(out$type_m_mcse))

  # With no true effect, all three are NA, as Type S and Type M are.
  spec$fixed$coefficients$effect <- 0
  null <- power_mixed(spec, n_sims = 3)
  expect_true(all(is.na(c(null$type_s_lo, null$type_s_hi, null$type_m_mcse))))
})
