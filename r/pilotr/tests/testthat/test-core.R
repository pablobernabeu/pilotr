# Core public-API tests: spec building, simulation, the JSON round-trip, power, and the
# response families. Kept fast (small n_sims, no lme4) so they run well within CRAN limits.

gaussian_between <- function(seed = 2024) {
  build_spec(list(
    name = "t", seed = seed, n_subject = 64, design_kind = "between",
    factor_name = "group", lev1 = "control", lev2 = "treatment",
    intercept = 100, effect = 5, family = "gaussian", resp_name = "", sigma = 10))
}

test_that("build_spec and simulate_design produce a well-formed between-subjects design", {
  spec <- gaussian_between()
  d <- simulate_design(spec)
  expect_s3_class(d, "data.frame")
  expect_equal(nrow(d), 64)
  expect_true(all(c("subject", "group", "score") %in% names(d)))
  expect_null(spec$random)
  expect_setequal(unique(d$group), c("control", "treatment"))
})

test_that("the same specification and seed reproduce identical data", {
  a <- simulate_design(gaussian_between(seed = 7))
  b <- simulate_design(gaussian_between(seed = 7))
  expect_identical(a, b)
})

test_that("spec_json is valid JSON and round-trips through load_spec", {
  spec <- gaussian_between()
  js <- spec_json(spec)
  expect_true(jsonlite::validate(js))
  tmp <- tempfile(fileext = ".json"); on.exit(unlink(tmp))
  writeLines(js, tmp)
  expect_equal(simulate_design(load_spec(tmp))$score, simulate_design(spec)$score)
})

test_that("power_design returns Type S / Type M and a plausible power", {
  r <- power_design(gaussian_between(), n_sims = 100)
  expect_true(all(c("n_sims", "power", "type_s", "type_m", "true_effect", "mean_estimate")
                  %in% names(r)))
  expect_gte(r$power, 0); expect_lte(r$power, 1)
  expect_gt(r$power, 0.2); expect_lt(r$power, 0.8)   # d = 0.5, n = 32/group
  expect_equal(r$true_effect, 5)
})

# design_conditions() recommends a condition in which every effect is zero, so this input is not
# hypothetical. Dividing by it used to report type_m = Inf and a type_s that had degenerated into
# the proportion of positive estimates.
test_that("a zero true effect leaves Type S and Type M undefined rather than infinite", {
  spec <- gaussian_between()
  spec$fixed$coefficients[[1]] <- 0
  r <- power_design(spec, n_sims = 100)
  expect_equal(r$true_effect, 0)
  expect_true(is.nan(r$type_s))
  expect_true(is.nan(r$type_m))
  expect_false(is.infinite(r$type_m))
  expect_gte(r$power, 0); expect_lte(r$power, 1)   # power itself is still reported
})

test_that("power_mixed also refuses to divide by a zero true effect", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  spec <- build_spec(list(name = "null", seed = 1, design_kind = "within",
                          include_items = FALSE, n_subject = 40,
                          factor_name = "cond", lev1 = "a", lev2 = "b",
                          intercept = 6, effect = 0, subj_int_sd = 0.12,
                          subj_slope_sd = 0, subj_corr = 0,
                          family = "gaussian", resp_name = "", sigma = 0.3))
  out <- power_mixed(spec, n_sims = 4)
  expect_true(out$n_returned > 0L)          # so the NAs below are the guard, not a failed fit
  expect_true(all(is.na(out$type_s)))
  expect_true(all(is.na(out$type_m)))
  expect_false(any(is.infinite(out$type_m)))
})

test_that("default_response_name covers every family", {
  expect_equal(default_response_name("gaussian"), "score")
  expect_equal(default_response_name("lognormal"), "RT")
  expect_equal(default_response_name("shifted_lognormal"), "RT")
  expect_equal(default_response_name("bernoulli"), "accuracy")
  expect_equal(default_response_name("poisson"), "count")
  expect_equal(default_response_name("ordinal"), "rating")
  expect_equal(default_response_name("beta"), "proportion")
})

test_that("each response family simulates a column on the expected scale", {
  fam <- function(family, extra = list()) {
    p <- c(list(name = "f", seed = 1, n_subject = 80, design_kind = "between",
                factor_name = "group", lev1 = "a", lev2 = "b", intercept = 0,
                effect = 0.5, family = family, resp_name = ""), extra)
    simulate_design(build_spec(p))
  }
  expect_true(all(fam("bernoulli")$accuracy %in% c(0, 1)))
  pois <- fam("poisson", list(intercept = 1.5, effect = 0.3))$count
  expect_true(all(pois >= 0) && all(pois == round(pois)))
  ord <- fam("ordinal", list(intercept = 0, effect = 0.8,
                             thresholds = "-2, -0.6, 0.6, 2"))$rating
  expect_true(all(ord >= 1) && all(ord <= 5))
  prop <- fam("beta", list(intercept = 0, effect = 0.8, phi = 8))$proportion
  expect_true(all(prop > 0) && all(prop < 1))
})

# The refusal is byte-identical to the Python twin's, which normalises its "inf" to R's
# "Inf" when the linear predictor overflows exp().
poisson_refusal <- function(mean) sprintf(paste0(
  "the poisson mean exp(eta) = %s is too large for the inverse-CDF sampler: beyond a mean of ",
  "about 708.4, exp(-mean) is subnormal and the walk loses accuracy. Lower the poisson ",
  "intercept or coefficients until the implied mean is at most 708."), mean)

test_that("a poisson mean beyond the sampler's accurate range is refused", {
  pois <- function(intercept, n = 4) build_spec(list(
    name = "p", seed = 1, n_subject = n, design_kind = "between",
    factor_name = "group", lev1 = "a", lev2 = "b", intercept = intercept,
    effect = 0, family = "poisson", resp_name = ""))
  # From a mean of about 738, the walk's total mass departs from 1. Over 40,000 draws, a mean of
  # 742 was refused for 4% of them. A mean of 744 drew counts with a mean of 733.5 and a
  # variance of 407.
  expect_error(simulate_design(pois(log(742))), poisson_refusal("742"), fixed = TRUE)
  # An intercept of 7 implies a mean of exp(7), about 1097, where exp(-mean) underflows to
  # zero; every count then used to come back as the sampler's 1e6 iteration cap.
  expect_error(simulate_design(pois(7)), poisson_refusal("1096.63"), fixed = TRUE)
  expect_error(simulate_design(pois(800)), poisson_refusal("Inf"), fixed = TRUE)
  # A mean of 700 is served, and its counts centre on it.
  d <- simulate_design(pois(log(700), n = 200))
  expect_true(all(d$count == round(d$count)))
  expect_lt(abs(mean(d$count) - 700), 5 * sqrt(700 / 200))
  # Feasible means are untouched: the shipped poisson example still draws real counts.
  d <- simulate_design(load_spec(pilotr_example("poisson_counts_between")))
  expect_true(all(d$count >= 0) && all(d$count < 1e6))
})

test_that("the poisson cap is the double nearest -log of the smallest normal double", {
  # 1022 * log(2) rounds to the literal both twins compare against. exp(-literal) is still a
  # normal double, so that mean is served. The next double up gives a subnormal and is refused.
  expect_identical(pilotr:::.poisson_inv(708.3964185322641, 0.5), 708)
  expect_error(pilotr:::.poisson_inv(0x1.6232bdd7abcd3p+9, 0.5), poisson_refusal("708.396"),
               fixed = TRUE)
})

test_that("the poisson walk evaluates each term as (p * lambda) / k, as the Python twin does", {
  # From exp(-3.7), the cumulative sum at k = 3 is 0.4941532441504185 with this bracketing and
  # 0.49415324415041845 with p * (lambda / k), which core.py used. This uniform therefore drew 3
  # in R and 4 in Python. The case was derived from the correctly rounded exp(-3.7).
  skip_if_not(identical(exp(-3.7), 0x1.9511fc6871044p-6),
              "this platform's exp(-3.7) is not the correctly rounded value")
  expect_identical(pilotr:::.poisson_inv(3.7, 0.4941532441504185), 3)
})

test_that("build_spec carries sigma through for the lognormal family", {
  spec <- build_spec(list(name = "ln", seed = 1, design_kind = "between",
                          factor_name = "g", lev1 = "a", lev2 = "b", n_subject = 10,
                          intercept = 0, effect = 0.5, family = "lognormal",
                          resp_name = "", sigma = 1))
  expect_equal(spec$response$sigma, 1)
  expect_equal(spec$response$name, "RT")
  d <- simulate_design(spec)
  expect_equal(nrow(d), 10)
  expect_true(all(d$RT > 0))
})

test_that("model_data() analyses both lognormal families on the log scale", {
  # The two families draw the same eta + sigma * z, so the analysis response is one column
  # whichever is declared. Python's power_mixed() follows the same rule, and it fitted a
  # lognormal response on its raw scale until it did.
  spec <- build_spec(list(name = "ln", seed = 1, design_kind = "within",
                          include_items = TRUE, n_subject = 12, n_item = 8,
                          factor_name = "cond", lev1 = "a", lev2 = "b",
                          intercept = 6, effect = 0.1, subj_int_sd = 0.12,
                          subj_slope_sd = 0, item_int_sd = 0.08, item_slope_sd = 0,
                          family = "lognormal", resp_name = "", sigma = 0.3))
  shifted <- spec
  shifted$response$family <- "shifted_lognormal"
  shifted$response$shift <- 200
  d <- simulate_design(spec)
  y <- model_data(spec, d)$.y
  expect_identical(y, log(d$RT))
  expect_equal(model_data(shifted, simulate_design(shifted))$.y, y, tolerance = 1e-9)
})

test_that("per_subject must lie between 1 and the number of items", {
  spec <- build_spec(list(name = "pc", seed = 1, design_kind = "within",
                          include_items = TRUE, n_subject = 2, n_item = 3,
                          factor_name = "cond", lev1 = "a", lev2 = "b",
                          intercept = 6, effect = 0.05, subj_int_sd = 0.1,
                          subj_slope_sd = 0, item_int_sd = 0.1, item_slope_sd = 0,
                          family = "gaussian", resp_name = "", sigma = 0.3))
  spec$units$item$per_subject <- 5
  expect_error(simulate_design(spec), "cannot exceed the number of items")
  spec$units$item$per_subject <- 0
  expect_error(simulate_design(spec), "at least 1")
})

# Until 0.3 power_mixed() refused any design without an item unit, and any design with more than one
# within-unit factor, because it fitted a formula written into the source rather than the one the
# specification implies. Both restrictions are gone, so what is checked now is that such a design
# runs and reports against the right effects.
test_that("power_mixed accepts a design with no item unit", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  spec <- build_spec(list(name = "w", seed = 1, design_kind = "within",
                          include_items = FALSE, n_subject = 40,
                          factor_name = "cond", lev1 = "a", lev2 = "b",
                          intercept = 6, effect = 0.4, subj_int_sd = 0.12,
                          subj_slope_sd = 0, subj_corr = 0,
                          family = "gaussian", resp_name = "", sigma = 0.3))
  out <- power_mixed(spec, n_sims = 4)
  expect_s3_class(out, "pilotr_power")
  expect_identical(names(out$power), "effect")
  expect_true(out$n_returned > 0L)
})

test_that("power_mixed accepts more than one within-unit factor", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  spec <- build_spec(list(name = "two", seed = 1, design_kind = "within",
                          include_items = TRUE, n_subject = 12, n_item = 8,
                          factor_name = "cond", lev1 = "a", lev2 = "b",
                          intercept = 6, effect = 0.05, subj_int_sd = 0.12,
                          subj_slope_sd = 0, item_int_sd = 0.08, item_slope_sd = 0,
                          family = "gaussian", resp_name = "", sigma = 0.3))
  spec$factors[[2]] <- list(name = "block", levels = c("x", "y"),
                            contrasts = list(blk = c(-0.5, 0.5)),
                            vary_within = c("subject", "item"))
  spec$fixed$coefficients$blk <- 0.03
  out <- power_mixed(spec, n_sims = 4)
  expect_setequal(names(out$power), c("effect", "blk"))
})

test_that("every reported rate carries its Monte Carlo error", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  spec <- build_spec(list(name = "mc", seed = 1, design_kind = "within",
                          include_items = TRUE, n_subject = 10, n_item = 8,
                          factor_name = "cond", lev1 = "a", lev2 = "b",
                          intercept = 6, effect = 0.05, subj_int_sd = 0.12,
                          subj_slope_sd = 0, item_int_sd = 0.08, item_slope_sd = 0,
                          family = "gaussian", resp_name = "", sigma = 0.3))
  out <- power_mixed(spec, n_sims = 6)
  expect_equal(unname(out$power_mcse),
               unname(sqrt(out$power * (1 - out$power) / out$n_returned)), tolerance = 1e-12)
  expect_true(out$power_lo <= out$power && out$power <= out$power_hi)
  # The plain standard error is zero at the boundary, where the Wilson interval still has width.
  expect_identical(pilotr:::.wilson(0, 60)[1], 0)
  expect_true(pilotr:::.wilson(0, 60)[2] > 0)
  expect_identical(pilotr:::.wilson(1, 60)[2], 1)
  expect_true(pilotr:::.wilson(1, 60)[1] < 1)
})

test_that("the fit counts distinguish singular fits from clean convergence", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  # A maximal crossed structure at a small sample size is singular in most replicates, which is the
  # case the old n_converged reported as fully converged.
  spec <- build_spec(list(name = "sing", seed = 1, design_kind = "within",
                          include_items = TRUE, n_subject = 8, n_item = 6,
                          factor_name = "cond", lev1 = "a", lev2 = "b",
                          intercept = 6, effect = 0.05, subj_int_sd = 0.12,
                          subj_slope_sd = 0.04, subj_corr = 0.2,
                          item_int_sd = 0.08, item_slope_sd = 0.02, item_corr = -0.1,
                          family = "shifted_lognormal", resp_name = "", sigma = 0.3, shift = 200))
  out <- power_mixed(spec, n_sims = 6)
  expect_identical(out$n_attempted, 6L)
  expect_true(out$n_converged <= out$n_returned)
  expect_true(out$n_singular > 0L)
  # lme4 records its singular-fit notice among the optimiser's messages. Counted as a warning, it
  # made n_warning equal n_singular in every run, so a fit that failed to converge could not be told
  # apart from one that was merely singular. Three of these six fits are singular, as before the
  # fix, and the warning count no longer includes them.
  expect_lt(out$n_warning, out$n_singular)
  expect_identical(out$n_singular, 3L)
})

test_that("a singular fit is not also recorded as a fit with warnings", {
  skip_if_not_installed("lme4")
  # Every subject has the same mean, so the subject variance is estimated at the boundary.
  d <- data.frame(subject = factor(rep(1:10, each = 4)), x = rep(c(-0.5, 0.5), 20))
  d$y <- 1 + 0.3 * d$x + rep(c(0.1, -0.1, -0.2, 0.2), 10)
  f <- pilotr:::.fit_lmer(y ~ x + (1 | subject), d)
  expect_true(f$singular)
  expect_false(any(grepl("singular", f$messages, ignore.case = TRUE)))
  expect_false(f$converged)
})

test_that("an interaction random slope reaches the linear predictor", {
  # The confound-free test: the slope KEY stays present in both runs and only its standard
  # deviation changes, so the number of draws consumed is identical and any difference in the data
  # is attributable to the slope itself. Before 0.3 both runs were bit-identical.
  spec <- load_spec(pilotr_example("reading_time_continuous"))
  spec$spec_version <- "0.3"
  spec$response$round <- NULL
  spec$random$subject$correlations <- NULL
  spec$random$subject$slopes[["SyntaxPC:age"]] <- 0
  zero <- simulate_design(spec)
  spec$random$subject$slopes[["SyntaxPC:age"]] <- 0.4
  nonzero <- simulate_design(spec)
  expect_false(identical(zero$reading_time_per_word, nonzero$reading_time_per_word))
})

test_that("a mistyped coefficient key is refused rather than silently zeroed", {
  spec <- load_spec(pilotr_example("crossed_mixed_rt"))
  spec$fixed$coefficients <- list(cnod = 0.05)
  expect_error(validate_spec(spec), "neither a contrast column nor a predictor")
})

test_that("a non-positive-definite random-effect covariance is an error", {
  spec <- load_spec(pilotr_example("crossed_mixed_rt"))
  spec$random$subject$slopes[["extra"]] <- 0.1
  spec$fixed$coefficients[["extra"]] <- 0
  spec$predictors <- list(list(name = "extra", varies_by = "item", mean = 0, sd = 1))
  spec$random$subject$correlations <- list(`intercept,cond` = 0.95,
                                           `intercept,extra` = 0.95,
                                           `cond,extra` = -0.95)
  expect_error(simulate_design(spec), "not positive definite")
})

# One specification file is run through both engines, so the same mistake has to be reported the
# same way by both. Each message below is byte-identical to the Python twin's.
test_that("a whole version is read the same however it was written", {
  spec <- load_spec(pilotr_example("between_2group_gaussian"))
  newer <- "declares spec_version 1.0, which is newer"
  # R renders the JSON number 1.0 as "1", Python as "1.0". Neither engine may call it malformed.
  spec$spec_version <- 1.0
  expect_error(validate_spec(spec), newer, fixed = TRUE)
  spec$spec_version <- "1"
  expect_error(validate_spec(spec), newer, fixed = TRUE)
})

test_that("a non-object unit is reported rather than crashing", {
  spec <- load_spec(pilotr_example("between_2group_gaussian"))
  spec$units$subject <- 5
  expect_error(validate_spec(spec), "'units.subject' must be an object", fixed = TRUE)
  # An empty object is a missing n, not a wrong shape, which is what the twin says too.
  spec$units$subject <- list()
  expect_error(validate_spec(spec), "'units.subject.n' must be a whole number", fixed = TRUE)
})

test_that("a non-whole seed truncates, as int(abs(seed)) does in the twin", {
  # Reachable only through validate = FALSE, the fast path the replicate loops use. Rounding here
  # handed the two engines different data from one specification.
  expect_identical(make_rng(2.7)$uniform(), make_rng(2)$uniform())
  expect_identical(make_rng(3.5)$uniform(), make_rng(3)$uniform())
  expect_identical(make_rng(-2.7)$uniform(), make_rng(2)$uniform())
})

test_that("replicate seeds are distinct and reproducible", {
  s <- replicate_seeds(90210, 500)
  expect_length(unique(s), 500)
  expect_identical(s, replicate_seeds(90210, 500))
  expect_true(all(s >= 1))
  # Not the old consecutive rule.
  expect_false(identical(s[1:3], c(90210, 90211, 90212)))
})

test_that("spec_json round-trips a coefficient exactly", {
  spec <- load_spec(pilotr_example("crossed_mixed_rt"))
  spec$fixed$coefficients$cond <- 1 / 3
  f <- tempfile(fileext = ".json")
  writeLines(spec_json(spec), f)
  back <- load_spec(f)
  expect_identical(as.numeric(back$fixed$coefficients$cond), 1 / 3)
  # And it stays readable: a value typed as 0.3 is not written as 0.29999999999999999.
  spec$fixed$coefficients$cond <- 0.3
  expect_true(grepl('"cond": 0.3', spec_json(spec), fixed = TRUE))
})

test_that("spec_json writes numbers a build without long doubles reads back", {
  spec <- load_spec(pilotr_example("crossed_mixed_rt"))
  spec$fixed$coefficients$cond <- 1 / 3
  txt <- spec_json(spec)
  # Sixteen significant digits, not the seventeen jsonlite writes on its own. The shorter form
  # is reached only by formatting the double itself. Shortening jsonlite's output instead meant
  # reading it back with as.numeric(), which without long doubles lands one unit in the last
  # place low, and the specification then held a coefficient the user never set.
  expect_true(grepl('"cond": 0.3333333333333333', txt, fixed = TRUE))
  expect_false(grepl("0.33333333333333331", txt, fixed = TRUE))

  # A build without long-double arithmetic (no_long_double(), in helper-numbers.R) reads a
  # decimal exactly when its digits stay within 2^53 and its exponent within +/-22, which holds
  # for the numbers of this document at their shortest widths, 1/3's sixteen digits included. It
  # would not hold for seventeen digits, and it does not hold for every sixteen-digit number.
  # Blank the string literals first, so that a factor level such as "group 2" is not read as a
  # number, then check every number the document actually carries.
  outside_strings <- gsub('"[^"]*"', '""', txt)
  toks <- regmatches(outside_strings,
                     gregexpr("-?[0-9][0-9.eE+-]*", outside_strings))[[1]]
  expect_gt(length(toks), 10)
  for (tok in toks) expect_identical(no_long_double(tok), as.numeric(tok))
})

# jsonlite reads a number with the C library's strtod(), which rounds correctly, and so does
# Python's json. On a build with long doubles, which include R on x86 Windows and Linux,
# as.numeric() rounds a decimal twice, first to long double and then to double. spec_json()
# used to check its 15- and 16-digit forms with as.numeric(), and accepted these three, which
# load_spec() and Python then read one unit in the last place away.
test_that("spec_json writes numbers that load_spec reads back as the same double", {
  vals <- c(-0x1.0efdc5258p+3, 0x1.e2b06cc2p+2, 0x1.1fe271f9p+3)
  spec <- load_spec(pilotr_example("between_2group_gaussian"))
  spec$fixed$intercept <- vals[1]
  spec$fixed$coefficients$grp <- vals[2]
  spec$response$sigma <- vals[3]
  spec$response$round <- NULL
  f <- tempfile(fileext = ".json")
  on.exit(unlink(f))
  writeLines(spec_json(spec), f)
  back <- load_spec(f)
  expect_identical(c(back$fixed$intercept, back$fixed$coefficients$grp, back$response$sigma),
                   vals)
  # With the response unrounded, the saved specification simulated different data in every
  # cell from the specification it was saved from.
  expect_identical(simulate_design(back), simulate_design(spec))
})

test_that("spec_json writes every double so that jsonlite reads it back exactly", {
  # Two kinds of double, made without the random-number generator. The first has the 32-bit
  # resolution of runif() over (-10, 10). The second has a full 53-bit mantissa and spans 31
  # orders of magnitude, so that the exponent forms are written too.
  i <- seq_len(3000)
  coarse <- ((i * 2654435761) %% 2^32) / 2^32 * 20 - 10
  fine <- sin(i) * 10^(i %% 31 - 15)
  vals <- c(coarse, fine, 1 / 3, 0.1 + 0.2, 1e23, -2.5e-17, 1e-300, .Machine$double.xmax)
  back <- jsonlite::fromJSON(spec_json(list(values = vals)))$values
  expect_identical(as.numeric(back), vals)
})

# Evaluate an exported script as a user would, apart from its library(pilotr) line: the
# package is already loaded wherever these tests run.
run_exported <- function(script) {
  env <- new.env(parent = asNamespace("pilotr"))
  for (e in parse(text = script, keep.source = FALSE))
    if (!identical(e, quote(library(pilotr)))) eval(e, env)
  env
}

# Nine numbers that generate_r_script() used to write as R literals which a build without long
# doubles reads one unit in the last place away. The six fitted values the verifier found need
# seventeen digits, or sixteen whose digits exceed 2^53. The other three are the values that
# spec_json() wrote one unit away for jsonlite.
test_that("generate_r_script carries every number of the specification exactly", {
  spec <- build_spec(list(
    name = "exact", seed = 11, n_subject = 6, n_item = 4, design_kind = "within",
    include_items = TRUE, factor_name = "cond", lev1 = "a", lev2 = "b (2)",
    intercept = -0x1.0efdc5258p+3, effect = 0x1.e2b06cc2p+2,
    subj_int_sd = 0.9923772805192401, subj_slope_sd = 1.1015114923556097,
    subj_corr = -0.9531743346778423, item_int_sd = 1.9022786173852906,
    item_slope_sd = 0.09470803828730423, item_corr = -0.48358814964064667,
    family = "gaussian", resp_name = "", sigma = 0x1.1fe271f9p+3))
  spec$response$round <- NULL
  numbers <- function(s) c(
    s$fixed$intercept, s$fixed$coefficients$effect, s$response$sigma,
    s$random$subject$intercept_sd, s$random$subject$slopes$effect,
    s$random$subject$correlations[["intercept,effect"]], s$random$item$intercept_sd,
    s$random$item$slopes$effect, s$random$item$correlations[["intercept,effect"]])
  script <- generate_r_script(spec)
  run <- run_exported(script)
  expect_identical(numbers(run$spec), numbers(spec))
  expect_identical(run$data, simulate_design(spec))
  # No number is left for R's own parser to read inexactly on some build.
  for (tok in r_decimal_tokens(script))
    expect_identical(no_long_double(tok), as.numeric(tok), info = tok)
})

# R's parser reads a decimal as as.numeric() does. On a build with long doubles, that rounds it
# twice, to 64 bits and then to 53, and the second rounding can go the wrong way when the first
# lands on the midpoint between two doubles. Correctly rounded, "-8.46847779583186" is
# -0x1.0efdc52580001p+3, which a build without long doubles reads, while an x86 build reads its
# neighbour. A script written on one kind of build can be run on the other, so whichever build
# writes it, such a decimal never reaches the script as R source.
test_that("a number is written as R source only where every build reads it alike", {
  z <- -0x1.0efdc52580001p+3
  expect_false(pilotr:::.r_reads_back("-8.46847779583186", z))
  lit <- pilotr:::.num_literal(z)
  expect_false(identical(lit, "-8.46847779583186"))
  expect_identical(eval(parse(text = lit)), z)
  # A value typed with up to four decimal places stays a plain decimal. The last three lie only
  # 1.6 units of a 64-bit significand from a midpoint, the least that four places allow.
  for (v in c(10, 0.05, 0.95, 1.01, 0.01, 0.0001, 200.2119, 613.4101, -898.0274))
    expect_identical(pilotr:::.num_literal(v), sprintf("%.15g", v), info = format(v))
})

test_that("a label written to imitate the number tag survives spec_json", {
  spec <- load_spec(pilotr_example("between_2group_gaussian"))
  # Numbers reach the document as tagged strings and are unwrapped afterwards, so a label
  # spelled like a tagged number would come back as a number rather than as itself. The tag is
  # lengthened until the specification's own text does not contain it.
  spec$response$name <- "@pilotr-number@1.5@pilotr-number@"
  spec$fixed$coefficients$grp <- 1 / 3
  f <- tempfile(fileext = ".json")
  writeLines(spec_json(spec), f)
  back <- load_spec(f, validate = FALSE)
  expect_identical(back$response$name, spec$response$name)
  expect_identical(as.numeric(back$fixed$coefficients$grp), 1 / 3)
})

test_that("brms_bridge maps the beta family to brms's Beta()", {
  spec <- build_spec(list(name = "b", seed = 1, design_kind = "between",
                          factor_name = "g", lev1 = "a", lev2 = "b", n_subject = 10,
                          intercept = 0, effect = 0.8, family = "beta",
                          resp_name = "", phi = 8))
  bridge <- brms_bridge(spec)
  expect_equal(bridge$family, "Beta()")
})

test_that("a crossed within-design has by-subject and by-item random slopes", {
  spec <- build_spec(list(
    name = "c", seed = 1, n_subject = 20, n_item = 12, design_kind = "within",
    include_items = TRUE, factor_name = "cond", lev1 = "x", lev2 = "y",
    intercept = 6, effect = 0.05, subj_int_sd = 0.12, subj_slope_sd = 0.04, subj_corr = 0.2,
    item_int_sd = 0.08, item_slope_sd = 0.02, item_corr = -0.1,
    family = "shifted_lognormal", resp_name = "", sigma = 0.3, shift = 200))
  d <- simulate_design(spec)
  expect_equal(nrow(d), 20 * 12 * 2)
  expect_false(is.null(spec$random$subject$slopes))
  expect_false(is.null(spec$random$item$slopes))
  expect_true(all(d$RT > 200))
})

test_that("generate_r_script emits a self-contained, runnable script", {
  rs <- generate_r_script(gaussian_between())
  expect_true(grepl("library(pilotr)", rs, fixed = TRUE))
  expect_true(grepl("simulate_design", rs, fixed = TRUE))
  expect_true(grepl("spec <-", rs, fixed = TRUE))
})

test_that("pilotr_example lists the bundled specs and resolves each to a file", {
  names <- pilotr_example()
  expect_type(names, "character")
  expect_true(length(names) >= 1)
  expect_true("between_2group_gaussian" %in% names)
  for (nm in names) {
    path <- pilotr_example(nm)
    expect_true(file.exists(path))
    # Every shipped example loads and simulates without error.
    spec <- load_spec(path)
    d <- simulate_design(spec)
    expect_s3_class(d, "data.frame")
    expect_gt(nrow(d), 0L)
  }
  # The .json extension is optional.
  expect_identical(
    pilotr_example("between_2group_gaussian"),
    pilotr_example("between_2group_gaussian.json")
  )
})

test_that("pilotr_example rejects unknown or malformed names", {
  expect_error(pilotr_example("no_such_example"), "Unknown example")
  expect_error(pilotr_example(c("a", "b")), "single example name")
})

# validate_spec permits n = 1, and stats::var() of one value is NA. The total used to count that
# NA as zero and still call itself the sum of the parts.
test_that("total is the sum of the parts for a single-row design", {
  spec <- build_spec(list(name = "one", seed = 1, design_kind = "between",
                          factor_name = "g", lev1 = "a", lev2 = "b", n_subject = 1,
                          intercept = 0, effect = 0.5, family = "gaussian",
                          resp_name = "", sigma = 1))
  spec$random <- list(subject = list(intercept_sd = 0.2))
  v <- response_variance(spec)
  expect_false(is.na(v$fixed))
  expect_identical(v$fixed, 0)            # one row has no across-row variation
  expect_equal(v$total, sum(unlist(v[setdiff(names(v), "total")])))
  # And the total is now something calibrate_response() can solve against.
  expect_equal(response_variance(calibrate_response(spec, 2))$total, 2)
})

test_that("total is the sum of the parts for every shipped example", {
  for (nm in pilotr_example()) {
    v <- response_variance(pilotr_example(nm))
    expect_equal(v$total, sum(unlist(v[setdiff(names(v), "total")])),
                 info = nm)
  }
})
