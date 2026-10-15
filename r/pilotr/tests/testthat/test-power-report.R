# power_report(): a sample-size justification in Markdown, from a power result.
#
# The report is rendered from a record of counts and strings, and the Python twin renders the same
# record into the same bytes. fixtures/power-report-record.json is one such record, the case
# "mixed_full" of tools/parity/report_cases.json, and fixtures/power-report-expected.md is its
# text. python/tests/fixtures holds the same two files, and tools/parity/report_cross.py renders
# every case of report_cases.json in both twins and requires identical text.

fixture <- function(name) test_path("fixtures", name)

# A file's bytes as one UTF-8 string, whatever the locale. A checkout on Windows can give the
# file CRLF line endings, which are read back as the newlines the report writes.
read_utf8 <- function(path) {
  raw <- readBin(path, "raw", file.info(path)$size)
  out <- rawToChar(raw)
  Encoding(out) <- "UTF-8"
  gsub("\r\n", "\n", out, fixed = TRUE)
}

# A power_mixed() result built by hand, which needs neither lme4 nor any fitting. The report reads
# only these fields. The second effect returned an estimate in 198 of the 200 fits.
fake_mixed <- function(fitter = "lmerTest::lmer, Satterthwaite t") {
  out <- list(
    n_sims = 200, alpha = 0.05, fitter = fitter,
    n_attempted = 200L, n_returned = 200L, n_converged = 151L, n_singular = 46L, n_warning = 3L,
    power = c(effect = 141 / 200, age = 37 / 198),
    power_mcse = c(effect = 0.0322, age = 0.0277),
    power_lo = c(effect = 0.638, age = 0.139), power_hi = c(effect = 0.763, age = 0.247),
    n_significant = c(effect = 141L, age = 37L),
    true_effect = c(effect = 0.06, age = -0.02),
    mean_estimate = c(effect = 0.061, age = -0.021),
    type_s = c(effect = 0, age = 2 / 37),
    type_s_lo = c(effect = 0, age = 0.015), type_s_hi = c(effect = 0.0265, age = 0.177),
    type_m = c(effect = 1.2734, age = 2.1), type_m_mcse = c(effect = 0.0213, age = 0.19),
    rung_shares = c(maximal = 0.24, zcp = 0.115, intercepts = 0.645), fallback = "singular")
  class(out) <- c("pilotr_power", "list")
  out
}

# The result of the same design with its effects set to zero.
fake_null <- function() {
  out <- fake_mixed()
  out$n_returned <- 199L
  out$power <- c(effect = 11 / 199, age = 9 / 199)
  out$n_significant <- c(effect = 11L, age = 9L)
  out$true_effect <- c(effect = 0, age = 0)
  out$type_s[] <- NA_real_
  out$type_m[] <- NA_real_
  out$rung_shares <- NULL
  out$fallback <- NULL
  out
}

# A solved curve, in the shape target_n() returns.
fake_solved <- function() {
  list(value = 51.37236, lo = 40.12, hi = 67.4, level = 0.95, target = 0.8, se = 0.4,
       dispersion = 1.3071, x = "n_subject", y = "power", transform = "sqrt",
       intercept = -2, slope = 0.4, n_points = 7L, x_min = 8, x_max = 56,
       n = 52, n_lo = 41, n_hi = 68)
}

# The report's line for a rate, built here from the counts with the package's own helpers, so
# that the test holds the report to those helpers and not to a copy of its formatting code.
rate_line <- function(count, n) {
  r <- pilotr:::.rate_with_error(count, n, "p")
  sprintf("%.3f (Monte Carlo standard error %.3f; 95%% Wilson interval %.3f to %.3f)",
          r$p, r$p_mcse, r$p_lo, r$p_hi)
}

two_group_spec <- function(n = 64, effect = 5) {
  build_spec(list(name = "d", seed = 1, design_kind = "between", factor_name = "group",
                  lev1 = "a", lev2 = "b", n_subject = n, intercept = 100, effect = effect,
                  family = "gaussian", resp_name = "", sigma = 10))
}

test_that("a stored record renders as the expected text, byte for byte", {
  record <- jsonlite::fromJSON(fixture("power-report-record.json"), simplifyVector = FALSE)
  got <- pilotr:::.report_render(record)
  expect_identical(got, read_utf8(fixture("power-report-expected.md")))
})

test_that("a power_mixed() result becomes the record the parity case holds", {
  spec_record <- jsonlite::fromJSON(fixture("power-report-record.json"), simplifyVector = FALSE)
  rec <- pilotr:::.power_report_record(
    fake_mixed(), solved = fake_solved(), spec_file = NULL,
    justification = c(effect = spec_record$effects[[1]]$justification),
    null_result = fake_null(), rope = NULL)
  # Everything but the software line and the specification, which come from this session and
  # from a file the record names without having it.
  rec$reproducibility <- spec_record$reproducibility
  expect_identical(pilotr:::.report_render(rec), read_utf8(fixture("power-report-expected.md")))
})

test_that("each rate is computed in the report from its counts", {
  rep <- power_report(fake_mixed())
  expect_s3_class(rep, "pilotr_report")
  expect_true(grepl(paste0("- Power: ", rate_line(141, 200), ", from 141 of 200 replicates"),
                    rep, fixed = TRUE))
  # The second effect's power is over the 198 fits that returned an estimate for it.
  expect_true(grepl(paste0("- Power: ", rate_line(37, 198), ", from 37 of 198 replicates"),
                    rep, fixed = TRUE))
  # Type S: 2 of the 37 significant replicates had the wrong sign.
  ci <- pilotr:::.type_s_interval(2, 37, -0.02)
  expect_true(grepl(sprintf("- Type S: 0.054 (95%% Wilson interval %.3f to %.3f)", ci[1], ci[2]),
                    rep, fixed = TRUE))
  expect_true(grepl("- Type M: 2.100 (Monte Carlo standard error 0.190)", rep, fixed = TRUE))
})

test_that("the primary analysis names the fitter, alpha and the sidedness", {
  rep <- power_report(fake_mixed())
  expect_true(grepl("- Model and test: lmerTest::lmer, Satterthwaite t\n", rep, fixed = TRUE))
  expect_true(grepl("- Significance level: alpha = 0.05, two-sided\n", rep, fixed = TRUE))
  expect_true(grepl("- Simplification rule: maximal, then zcp, then intercepts, moving on from singular fits",
                    rep, fixed = TRUE))
  one <- power_report(fake_mixed("lme4::glmer (binomial), Wald z, one-sided (greater)"))
  expect_true(grepl("- Significance level: alpha = 0.05, one-sided (greater)\n", one, fixed = TRUE))
  less <- power_report(fake_mixed("stats::lm, t, one-sided (less)"))
  expect_true(grepl("alpha = 0.05, one-sided (less)\n", less, fixed = TRUE))
})

test_that("a power_design() result is reported over all its replicates", {
  res <- power_design(two_group_spec(), n_sims = 40, alternative = "greater")
  rep <- power_report(res, justification = "The smallest difference worth the cost of the course.")
  expect_true(grepl("### The difference between the two groups", rep, fixed = TRUE))
  expect_true(grepl(sprintf("- Power: %.3f (Monte Carlo standard error %.3f; 95%% Wilson interval %.3f to %.3f), from %d of 40 replicates",
                            res$power, res$power_mcse, res$power_lo, res$power_hi,
                            res$n_significant), rep, fixed = TRUE))
  expect_true(grepl("- Model and test: stats::t.test, two-sample t with pooled variance", rep,
                    fixed = TRUE))
  expect_true(grepl("alpha = 0.05, one-sided (greater)", rep, fixed = TRUE))
  expect_true(grepl("- Replicates: 40 simulated\n", rep, fixed = TRUE))
  expect_true(grepl("- Justification: The smallest difference worth the cost of the course.",
                    rep, fixed = TRUE))
  expect_false(grepl("- Fits:", rep, fixed = TRUE))
})

test_that("a solved curve gives the sample size with its interval, transform and dispersion", {
  rep <- power_report(fake_mixed(), solved = fake_solved())
  expect_true(grepl("## Desired power and sample size", rep, fixed = TRUE))
  expect_true(grepl("- Desired power: 0.8\n", rep, fixed = TRUE))
  expect_true(grepl("- Sample size: 52 (95% confidence interval 41 to 68), the exact solution 51.37 rounded up",
                    rep, fixed = TRUE))
  expect_true(grepl("the square root of `n_subject`", rep, fixed = TRUE))
  expect_true(grepl("- Heterogeneity factor: 1.31", rep, fixed = TRUE))

  # A real solve, which target_n() returns with its whole-number fields.
  curve <- data.frame(n_subject = c(20, 40, 60, 80, 120, 160),
                      power = c(0.21, 0.38, 0.52, 0.63, 0.78, 0.87),
                      n_returned = c(200, 200, 190, 200, 180, 200))
  solved <- target_n(curve, target = 0.8)
  rep <- power_report(fake_mixed(), solved = solved)
  expect_true(grepl(sprintf("- Sample size: %d (95%% confidence interval %d to %d), the exact solution %.2f rounded up",
                            solved$n, solved$n_lo, solved$n_hi, solved$value), rep, fixed = TRUE))

  # solve_curve() over an effect size has no whole-number fields.
  effects <- data.frame(cond = c(0.1, 0.2, 0.3, 0.4, 0.5), power = c(0.11, 0.31, 0.58, 0.79, 0.92),
                        n_sims = 400)
  mde <- solve_curve(effects, target = 0.8, transform = "identity")
  rep <- power_report(fake_mixed(), solved = mde)
  expect_true(grepl("## Desired power and the solved value", rep, fixed = TRUE))
  expect_true(grepl(sprintf("- Solved value of `cond`: %g (95%% confidence interval %g to %g)",
                            mde$value, mde$lo, mde$hi), rep, fixed = TRUE))
  expect_error(power_report(fake_mixed(), solved = list(value = 3)),
               "`solved` must be the result of target_n() or solve_curve().", fixed = TRUE)
})

test_that("the null result gives the rejection rate when every focal effect is zero", {
  rep <- power_report(fake_mixed(), null_result = fake_null())
  expect_true(grepl("## Power when the null hypothesis is true", rep, fixed = TRUE))
  expect_true(grepl(paste0("- `effect`: ", rate_line(11, 199), ", from 11 of 199 replicates"),
                    rep, fixed = TRUE))
  expect_true(grepl("- Null condition: 200 simulated, 199 returned an estimate", rep, fixed = TRUE))
  # A result whose effects are not zero is not a null result.
  expect_error(power_report(fake_mixed(), null_result = fake_mixed()),
               "`null_result` must come from a design in which every focal effect is zero. Its effect 'effect' has a true value of 0.06.",
               fixed = TRUE)
})

test_that("justifications attach to their effects, and a missing one is said to be missing", {
  rep <- power_report(fake_mixed(), justification = c(age = "Ageing studies report this."))
  expect_true(grepl("### `age`\n\n- True value: -0.02\n- Justification: Ageing studies report this.",
                    rep, fixed = TRUE))
  expect_true(grepl("### `effect`\n\n- True value: 0.06\n- Justification: not given", rep,
                    fixed = TRUE))
  both <- power_report(fake_mixed(), justification = "One reason for both.")
  expect_identical(lengths(regmatches(both, gregexpr("One reason for both.", both, fixed = TRUE))),
                   2L)
  expect_error(power_report(fake_mixed(), justification = c(cond = "x")),
               "`justification` names 'cond', which is not a focal effect of `x`. Its focal effects are 'effect' and 'age'.",
               fixed = TRUE)
  form <- "`justification` must be one string, which applies to every focal effect, or strings named by focal effect."
  expect_error(power_report(fake_mixed(), justification = c("a", "b")), form, fixed = TRUE)
  expect_error(power_report(fake_mixed(), justification = 1), form, fixed = TRUE)
  res <- power_design(two_group_spec(), n_sims = 10)
  expect_error(power_report(res, justification = c(group = "x")),
               "`x` has one effect, without a name, so `justification` must be a single string.",
               fixed = TRUE)
})

test_that("the specification file gives the checksum of its bytes, its seed and its units", {
  spec <- two_group_spec()
  path <- tempfile(fileext = ".json")
  on.exit(unlink(path), add = TRUE)
  writeLines(spec_json(spec), path)
  res <- power_design(spec, n_sims = 10)
  rep <- power_report(res, spec_file = path)
  expect_true(grepl(sprintf("- Specification: %s, MD5 checksum %s\n", basename(path),
                            unname(tools::md5sum(path))), rep, fixed = TRUE))
  expect_true(grepl("- Seed: 1\n", rep, fixed = TRUE))
  expect_true(grepl("- Simulated units: subject 64\n", rep, fixed = TRUE))
  expect_true(grepl(sprintf("- Software: pilotr %s, R %s\n", utils::packageVersion("pilotr"),
                            getRversion()), rep, fixed = TRUE))
  none <- power_report(res)
  expect_true(grepl("- Specification: not given, so neither its seed nor its checksum is recorded",
                    none, fixed = TRUE))
  missing <- file.path(tempdir(), "no-such-spec.json")
  expect_error(power_report(res, spec_file = missing),
               sprintf("`spec_file` must be the path of the specification file the analysis read. No file is found at '%s'.",
                       missing), fixed = TRUE)
})

test_that("the software line names lme4 and lmerTest when the fitter used them", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  rep <- power_report(fake_mixed())
  expect_true(grepl(sprintf("- Software: pilotr %s, R %s, lme4 %s, lmerTest %s\n",
                            utils::packageVersion("pilotr"), getRversion(),
                            utils::packageVersion("lme4"), utils::packageVersion("lmerTest")),
                    rep, fixed = TRUE))
})

test_that("a curve, or anything else that is not a result at one design, is refused", {
  curve <- data.frame(n_subject = c(10, 20), effect = "effect", power = c(0.2, 0.4),
                      n_returned = 50)
  msg <- "`x` must be the result of power_mixed(), power_design() or precision_design() at one design. A curve goes to target_n(), and its result to `solved`."
  expect_error(power_report(curve), msg, fixed = TRUE)
  expect_error(power_report(list(a = 1)), msg, fixed = TRUE)
  expect_error(power_report(fake_mixed(), rope = 0.1),
               "`rope` states the region of practical equivalence of a precision_design() result.",
               fixed = TRUE)
})

test_that("an effect without a true value or without estimates is reported as such", {
  x <- fake_mixed()
  x$true_effect[] <- NA_real_
  x$type_s[] <- NA_real_
  x$type_m[] <- NA_real_
  x$type_m_mcse[] <- NA_real_
  x$power[["age"]] <- NA_real_
  x$n_significant[["age"]] <- 0L
  rep <- power_report(x)
  expect_true(grepl("- True value: not given, so Type S and Type M are not defined", rep,
                    fixed = TRUE))
  expect_true(grepl("- Power: not available, since no replicate returned an estimate for this effect",
                    rep, fixed = TRUE))
  # Type M with a single significant replicate has no standard error.
  y <- fake_mixed()
  y$type_m_mcse[["age"]] <- NA_real_
  expect_true(grepl("- Type M: 2.100 (too few significant replicates for a Monte Carlo standard error)",
                    power_report(y), fixed = TRUE))
})

test_that("a joint test lists each set's coefficients and has no Type S or Type M", {
  x <- fake_mixed("lmerTest::lmer, Satterthwaite F")
  keep <- "effect"
  for (k in c("power", "power_mcse", "power_lo", "power_hi", "n_significant", "true_effect",
              "mean_estimate", "type_s", "type_s_lo", "type_s_hi", "type_m", "type_m_mcse"))
    x[[k]] <- stats::setNames(x[[k]][keep], "interaction")
  x$true_effect[] <- NA_real_
  x$type_s[] <- NA_real_
  x$type_m[] <- NA_real_
  x$sets <- list(interaction = c("A1_B", "A2_B"))
  rep <- power_report(x)
  expect_true(grepl("- Coefficients tested jointly: A1_B, A2_B", rep, fixed = TRUE))
  expect_true(grepl("alpha = 0.05, a joint test of each set", rep, fixed = TRUE))
  expect_false(grepl("Type S", rep, fixed = TRUE))
  expect_false(grepl("True value", rep, fixed = TRUE))
})

test_that("a precision_design() result gives the interval and ROPE decisions", {
  prec <- data.frame(param = c("effect", "age"), true = c(0.06, -0.02),
                     mean_ci_width = c(0.0812, 0.0411),
                     p_meaningful = c(9 / 20, 2 / 20), p_meaningful_mcse = 0,
                     p_meaningful_lo = 0, p_meaningful_hi = 0,
                     p_equivalent = c(0, 3 / 20), p_equivalent_mcse = 0,
                     p_equivalent_lo = 0, p_equivalent_hi = 0,
                     n_attempted = 20L, n_returned = 20L, n_converged = 18L, n_singular = 2L,
                     n_warning = 0L, fitter = "lme4::lmer, Wald z", stringsAsFactors = FALSE)
  rep <- power_report(prec, rope = 0.02, justification = c(effect = "A pilot estimate."))
  expect_true(grepl("- Model and interval: lme4::lmer, Wald z", rep, fixed = TRUE))
  expect_true(grepl("a 95% Wald interval wholly outside the region of practical equivalence, from -0.02 to 0.02",
                    rep, fixed = TRUE))
  expect_true(grepl(paste0("- Interval wholly outside the region (a meaningful effect): ",
                           rate_line(9, 20), ", from 9 of 20 replicates"), rep, fixed = TRUE))
  expect_true(grepl(paste0("- Interval wholly inside the region (practical equivalence): ",
                           rate_line(3, 20), ", from 3 of 20 replicates"), rep, fixed = TRUE))
  expect_true(grepl("- Mean interval width: 0.0812", rep, fixed = TRUE))
  expect_false(grepl("Significance level", rep, fixed = TRUE))
  # Without `rope`, the report says where the region was set.
  unstated <- power_report(prec)
  expect_true(grepl("the region of practical equivalence given to precision_design()", unstated,
                    fixed = TRUE))
  # A target_n() result over p_meaningful is a desired probability, not a desired power.
  solved <- fake_solved()
  solved$y <- "p_meaningful"
  solved$target <- 0.9
  rep <- power_report(prec, solved = solved, rope = c(-0.01, 0.03))
  expect_true(grepl("- Desired probability of a meaningful decision: 0.9", rep, fixed = TRUE))
  expect_true(grepl("from -0.01 to 0.03", rep, fixed = TRUE))
  expect_error(power_report(prec, null_result = fake_null()),
               "`null_result` must be a result of the same kind as `x`.", fixed = TRUE)
  expect_error(power_report(prec, rope = c(0.1, -0.1)), "`rope` must be", fixed = TRUE)
})

test_that("printing a report writes its Markdown and nothing else", {
  rep <- power_report(fake_mixed())
  txt <- as.character(unclass(rep))
  expect_identical(substr(txt, 1, 28), "# Sample-size justification\n")
  expect_true(endsWith(txt, "\n") && !endsWith(txt, "\n\n"))
  out <- utils::capture.output(shown <- withVisible(print(rep)))
  expect_identical(paste0(paste(out, collapse = "\n"), "\n"), txt)
  expect_false(shown$visible)
  expect_identical(shown$value, rep)
})

test_that("the report is the same on every call", {
  expect_identical(power_report(fake_mixed(), solved = fake_solved(), null_result = fake_null()),
                   power_report(fake_mixed(), solved = fake_solved(), null_result = fake_null()))
})
