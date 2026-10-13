# A print method has to reach the reader as a single block. cli and message()
# write to the message stream, which knitr collects separately from standard
# output, so a method that put its header on one stream and its table on the
# other rendered one printed object as two boxes on the documentation site.
# print.pilotr_power writes throughout with cat() and print(); these tests hold
# it there. print.pilotr_bridge, the package's other method, is tested in
# test-brms-bridge.R alongside the function whose result it displays.

# A power result built by hand rather than by power_mixed(), which needs lme4
# and a few seconds of fitting. The print method only reads these fields.
fake_power <- function() {
  out <- list(
    n_sims = 200L, alpha = 0.05,
    n_attempted = 200L, n_returned = 198L, n_converged = 190L,
    n_singular = 6L, n_warning = 2L,
    power = c(effect = 0.62), power_mcse = c(effect = 0.034),
    power_lo = c(effect = 0.551), power_hi = c(effect = 0.686),
    n_significant = c(effect = 123L), true_effect = c(effect = 0.4),
    mean_estimate = c(effect = 0.41),
    type_s = c(effect = 0), type_s_lo = c(effect = 0), type_s_hi = c(effect = 0.0303),
    type_m = c(effect = 1.32), type_m_mcse = c(effect = 0.0214))
  class(out) <- c("pilotr_power", "list")
  out
}

test_that("printing a power result writes nothing to the message stream", {
  # Sink standard output so the table itself does not clutter the reporter.
  con <- file(nullfile(), open = "wt")
  sink(con)
  on.exit({ sink(); close(con) }, add = TRUE)
  expect_length(capture.output(print(fake_power()), type = "message"), 0L)
})

test_that("the header and the table share standard output", {
  out <- capture.output(print(fake_power()))
  expect_true(any(grepl("Simulation-based power over 200 replicates", out)))
  expect_true(any(grepl("198 returned", out)))
  # The table has to be in the same capture as the header, not a second stream.
  expect_true(any(grepl("^\\s*effect\\s", out)))
  expect_true(any(grepl("type_m", out)))
})

test_that("the singularity note is part of the same block", {
  x <- fake_power()
  x$n_singular <- 120L
  out <- capture.output(print(x))
  note <- grep("boundary-singular", out)
  expect_length(note, 1L)
  expect_true(note < length(out))
})

test_that("the fitter and its test are printed with the fit counts", {
  x <- fake_power()
  x$fitter <- "lme4::glmer (binomial), Wald z"
  out <- capture.output(print(x))
  line <- grep("lme4::glmer (binomial), Wald z", out, fixed = TRUE)
  expect_length(line, 1L)
  expect_true(line < grep("198 returned", out))
  # A result saved by an earlier version has no fitter, and prints as it did.
  expect_false(any(grepl("fitter", capture.output(print(fake_power())))))
})

test_that("Type S and Type M are printed with their Monte Carlo error", {
  out <- capture.output(print(fake_power()))
  heading <- grep("Type S and Type M, over the significant replicates", out, fixed = TRUE)
  expect_length(heading, 1L)
  expect_true(heading > grep("[0.551, 0.686]", out, fixed = TRUE))
  expect_true(any(grepl("^\\s*effect\\s+type_s\\s+type_s_ci95\\s+type_m\\s+type_m_mcse\\s*$", out)))
  expect_true(any(grepl("[0.000, 0.030]", out, fixed = TRUE)))
  expect_true(any(grepl("0.021", out, fixed = TRUE)))
})

test_that("digits sets the decimal places of the rates and of the intervals", {
  out <- capture.output(print(fake_power(), digits = 2))
  expect_true(any(grepl("[0.55, 0.69]", out, fixed = TRUE)))
  expect_true(any(grepl("[0.00, 0.03]", out, fixed = TRUE)))
  expect_true(any(grepl(" 0.62 ", out, fixed = TRUE)))
  expect_false(any(grepl("0.551", out, fixed = TRUE)))
})

test_that("a result saved before the Type S interval existed prints without it", {
  x <- fake_power()
  x$type_s_lo <- x$type_s_hi <- x$type_m_mcse <- NULL
  out <- capture.output(print(x))
  expect_true(any(grepl("^\\s*effect\\s+type_s\\s+type_m\\s*$", out)))
  expect_false(any(grepl("type_s_ci95|type_m_mcse", out)))
})

test_that("printing returns its input invisibly", {
  x <- fake_power()
  con <- file(nullfile(), open = "wt")
  sink(con)
  on.exit({ sink(); close(con) }, add = TRUE)
  expect_invisible(print(x))
})
