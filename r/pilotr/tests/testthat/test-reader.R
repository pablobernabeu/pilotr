# How a specification is read: every field by its exact name, repeated keys refused, seeds kept
# within the range both twins read exactly, and malformed shapes reported in the validator's
# words. The Python twin's tests/test_reader.py covers the same ground from its side.

reader_between <- function(seed = 2024) {
  build_spec(list(
    name = "t", seed = seed, n_subject = 24, design_kind = "between",
    factor_name = "group", lev1 = "control", lev2 = "treatment",
    intercept = 100, effect = 5, family = "gaussian", resp_name = "", sigma = 10))
}

# A Gaussian design with one extra grouping factor and no `subject` or `item` random entry, the
# shape in which R's `$` used to read `random$subject` off an entry named, say, `subject_site`.
extra_group_spec <- function(gname, over = "subject") {
  list(
    name = "x", seed = 4417,
    units = list(subject = list(n = 12L), item = list(n = 4L)),
    factors = list(list(name = "condition", levels = c("a", "b"),
                        contrasts = list(cond = c(-0.5, 0.5)),
                        vary_within = c("subject", "item"))),
    fixed = list(intercept = 10, coefficients = list(cond = 0.5)),
    random = stats::setNames(list(list(over = over, n = 3L, intercept_sd = 2)), gname),
    response = list(family = "gaussian", name = "score", sigma = 1))
}

test_that("a grouping factor whose name extends 'subject' or 'item' adds no other effects", {
  a <- simulate_design(extra_group_spec("subject_site"))
  b <- simulate_design(extra_group_spec("site"))
  expect_identical(a$score, b$score)
  expect_identical(a$subject_site, b$site)

  a <- simulate_design(extra_group_spec("item_list", over = "item"))
  b <- simulate_design(extra_group_spec("cluster", over = "item"))
  expect_identical(a$score, b$score)
  expect_identical(a$item_list, b$cluster)
})

test_that("with strict = FALSE a misspelt field is ignored, as the warning says", {
  base <- reader_between()
  ref <- simulate_design(base)

  s <- base
  s$response$rounding <- 1
  expect_warning(validate_spec(s, strict = FALSE), "response.rounding", fixed = TRUE)
  expect_identical(simulate_design(s, validate = FALSE), ref)

  s <- base
  s$units$items <- list(n = 3L)
  expect_warning(validate_spec(s, strict = FALSE), "unknown unit 'items'", fixed = TRUE)
  expect_identical(simulate_design(s, validate = FALSE), ref)

  with_x <- base
  with_x$predictors <- list(list(name = "x", varies_by = "subject"))
  s <- with_x
  s$predictors[[1]]$mean_centre <- 50
  expect_warning(validate_spec(s, strict = FALSE), "predictors[1].mean_centre", fixed = TRUE)
  expect_identical(simulate_design(s, validate = FALSE), simulate_design(with_x))
})

test_that("a repeated key is refused, whatever `validate` says", {
  msg <- "the specification repeats the key 'effect' within one object; JSON leaves a repeated key undefined, and R and Python read it differently"
  js <- sub('"coefficients": {', '"coefficients": {"effect": 50, ', spec_json(reader_between()),
            fixed = TRUE)
  f <- tempfile(fileext = ".json"); on.exit(unlink(f))
  writeLines(js, f)
  expect_error(load_spec(f), msg, fixed = TRUE)
  expect_error(load_spec(f, validate = FALSE), msg, fixed = TRUE)
  expect_error(simulate_design(f, validate = FALSE), msg, fixed = TRUE)

  writeLines(sub("{", '{"seed": 7, ', spec_json(reader_between()), fixed = TRUE), f)
  expect_error(load_spec(f, validate = FALSE), "repeats the key 'seed'", fixed = TRUE)

  s <- reader_between()
  s$fixed$coefficients <- list(effect = 5, effect = 50)
  expect_error(validate_spec(s), msg, fixed = TRUE)
})

test_that("a seed must lie within the range both twins read exactly", {
  msg <- "'seed' must be a whole number between -9007199254740991 and 9007199254740991 (2^53 - 1), the range in which R and Python read a JSON integer exactly"
  s <- reader_between()
  s$seed <- 2^53
  expect_error(validate_spec(s), msg, fixed = TRUE)
  s$seed <- -2^53
  expect_error(validate_spec(s), msg, fixed = TRUE)
  s$seed <- 1.5
  expect_error(validate_spec(s), msg, fixed = TRUE)
  s$seed <- 2^53 - 1
  expect_silent(validate_spec(s))
  s$seed <- -(2^53 - 1)
  expect_silent(validate_spec(s))
})

test_that("build_spec() keeps a large seed and leaves a fractional one to the validator", {
  expect_identical(reader_between(seed = 2024)$seed, 2024L)
  big <- reader_between(seed = 3e9)
  expect_identical(big$seed, 3e9)
  expect_equal(nrow(simulate_design(big)), 24)
  frac <- reader_between(seed = 1.5)
  expect_identical(frac$seed, 1.5)
  expect_error(validate_spec(frac), "'seed' must be a whole number", fixed = TRUE)
})

test_that("a null level or vary_within entry is refused, as the Python twin refuses it", {
  s <- reader_between()
  s$factors[[1]]$levels <- c("control", NA)
  expect_error(validate_spec(s), "factors[1].levels must be an array of at least two strings",
               fixed = TRUE)

  s <- extra_group_spec("site")
  s$factors[[1]]$vary_within <- c("subject", NA)
  expect_error(validate_spec(s),
               "factors[1].vary_within must be a unit name or an array of unit names",
               fixed = TRUE)
})

test_that("malformed shapes are reported in the validator's own words", {
  s <- reader_between()
  s$response <- "gaussian"
  expect_error(validate_spec(s), "'response' must be an object", fixed = TRUE)

  s <- reader_between()
  s$predictors <- list("x")
  expect_error(validate_spec(s), "predictors[1] must be an object", fixed = TRUE)

  s <- reader_between()
  s$random <- list(subject = 5)
  expect_error(validate_spec(s), "random.subject must be an object", fixed = TRUE)

  s <- extra_group_spec("site")
  s$random$site$slopes <- 5
  expect_error(validate_spec(s), "random.site.slopes must be an object", fixed = TRUE)
})

test_that("a spec_version that is not one string or number is refused in the validator's words", {
  # An empty array or object stopped with "subscript out of bounds", and a two-element array was
  # read by its first element, where the Python twin called it malformed.
  msg <- "'spec_version' must be a single string of the form 'major.minor'"
  s <- reader_between()
  for (v in list(list(), character(0), c("0.3", "0.4"), NA, TRUE)) {
    s$spec_version <- v
    expect_error(validate_spec(s), msg, fixed = TRUE)
  }
})

test_that("a null correlated flag is refused and an array-valued varies_by reported once", {
  # `"correlated": [null]` reads as NA, which passed as a logical and meant FALSE.
  s <- extra_group_spec("site")
  s$spec_version <- "0.3"
  s$random$site$correlated <- NA
  expect_error(validate_spec(s), "random.site.correlated must be TRUE or FALSE", fixed = TRUE)

  # The value was pasted into the message, which a two-element array turned into two problems.
  s <- reader_between()
  s$predictors <- list(list(name = "x", varies_by = c("subject", "item")))
  s$fixed$coefficients$x <- 0.1
  err <- tryCatch(validate_spec(s), error = conditionMessage)
  expect_identical(
    lengths(regmatches(err, gregexpr("varies_by must be", err, fixed = TRUE))), 1L)
  expect_false(grepl(", not '", err, fixed = TRUE))
})

test_that("load_spec() reads a file with a byte-order mark and non-ASCII labels", {
  s <- reader_between()
  s$factors[[1]]$levels <- c("fácil", "Łatwy")
  f <- tempfile(fileext = ".json"); on.exit(unlink(f))
  con <- file(f, "wb")
  writeBin(as.raw(c(0xEF, 0xBB, 0xBF)), con)
  writeBin(charToRaw(enc2utf8(spec_json(s))), con)
  close(con)
  expect_silent(loaded <- load_spec(f))
  expect_identical(loaded$factors[[1]]$levels, c("fácil", "Łatwy"))
  expect_identical(simulate_design(loaded), simulate_design(s))
})

test_that("no specification field is read by partial matching", {
  old <- options(warnPartialMatchDollar = TRUE)
  on.exit(options(old))
  promote <- function(w) {
    if (grepl("partial match", conditionMessage(w), fixed = TRUE))
      stop("partial match: ", conditionMessage(w), call. = FALSE)
  }
  specs <- lapply(pilotr_example(), function(ex) load_spec(pilotr_example(ex)))
  specs <- c(specs, list(extra_group_spec("subject_site"),
                         extra_group_spec("item_list", over = "item")))
  for (spec in specs) {
    expect_no_error(withCallingHandlers({
      d <- simulate_design(spec)
      model_data(spec, d)
      model_formula(spec)
      brms_bridge(spec)
      response_variance(spec)
    }, warning = promote))
  }
})
