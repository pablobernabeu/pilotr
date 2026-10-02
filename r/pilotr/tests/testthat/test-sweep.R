# The tests in this file check what each grid point of sweep_spec() receives and what
# design_conditions() builds for it. Most of the analyses passed as `fn` only report the
# specification they were given, so no mixed model is fitted.

reading_spec <- function() load_spec(pilotr_example("reading_time_continuous"))

# One row holding the coefficients a grid point received, named as in the specification.
coefficient_row <- function(spec) as.data.frame(spec$fixed$coefficients, check.names = FALSE)

never <- function(spec) stop("`fn` was reached", call. = FALSE)

test_that("a design_conditions() sweep keeps every coefficient, in the specification's order", {
  spec <- reading_spec()
  coefs <- spec$fixed$coefficients
  seen <- sweep_spec(spec, "fixed$coefficients", design_conditions(cond = c(0.02, 0.04)),
                     coefficient_row)
  expect_identical(names(seen), c("coefficients", names(coefs)))
  expect_identical(seen$coefficients, 1:3)   # a list sweep records the grid index
  expect_identical(seen$cond, c(0, 0.02, 0.04))
})

test_that("the null condition zeroes only the effects the conditions name", {
  spec <- reading_spec()
  coefs <- spec$fixed$coefficients
  seen <- sweep_spec(spec, "fixed$coefficients", design_conditions(cond = c(0.02, 0.04)),
                     coefficient_row)
  others <- setdiff(names(coefs), "cond")
  for (k in others) expect_identical(seen[[k]], rep(coefs[[k]], 3), info = k)
  expect_identical(seen$cond[1], 0)
})

test_that("a .base value holds in every condition, the null one included", {
  spec <- reading_spec()
  coefs <- spec$fixed$coefficients
  seen <- sweep_spec(spec, "fixed$coefficients",
                     design_conditions(cond = 0.02, .base = list(age = 0.1)), coefficient_row)
  expect_identical(seen$age, c(0.1, 0.1))
  expect_identical(seen$cond, c(0, 0.02))
  expect_identical(seen$SyntaxPC, rep(coefs$SyntaxPC, 2))
})

test_that("a condition simulates exactly the data of the specification edited by hand", {
  spec <- reading_spec()
  # Rounding to five digits would hide the last-bit changes a different summation order makes.
  spec$response$round <- NULL
  by_hand <- spec
  by_hand$fixed$coefficients$cond <- 0.02
  want <- simulate_design(by_hand)$reading_time_per_word

  got <- sweep_spec(spec, "fixed$coefficients", design_conditions(cond = 0.02, .null = FALSE),
                    simulate_design)
  expect_identical(got$reading_time_per_word, want)

  # A .base listing the coefficients in another order leaves the specification's order in place.
  shuffled <- rev(spec$fixed$coefficients)
  got <- sweep_spec(spec, "fixed$coefficients",
                    design_conditions(cond = 0.02, .null = FALSE, .base = shuffled),
                    simulate_design)
  expect_identical(got$reading_time_per_word, want)
})

test_that("a condition naming a coefficient the specification lacks is refused", {
  spec <- reading_spec()
  expect_error(
    sweep_spec(spec, "fixed$coefficients", design_conditions(cnod = 0.02), never),
    paste("design_conditions() names 'cnod', which is not a coefficient of this specification;",
          "its coefficients are 'SyntaxPC', 'CoherencePC', 'age', 'cond', 'SyntaxPC:age',",
          "'CoherencePC:age', 'cond:age'"),
    fixed = TRUE)
  expect_error(
    sweep_spec(spec, "fixed$coefficients", design_conditions(cnod = 0.02, agee = 0.1), never),
    "design_conditions() names 'cnod', 'agee', which are not coefficients of this specification;",
    fixed = TRUE)

  none <- spec
  none$fixed$coefficients <- stats::setNames(list(), character(0))
  expect_error(
    sweep_spec(none, "fixed$coefficients", design_conditions(cond = 0.02), never),
    "its coefficients are (none)", fixed = TRUE)
})

test_that("design_conditions() refuses a path that is not a named list", {
  spec <- reading_spec()
  expect_error(
    sweep_spec(spec, "fixed$intercept", design_conditions(cond = 0.02), never),
    paste("design_conditions() merges its conditions into a named list such as",
          "fixed$coefficients, and `path` addresses fixed$intercept, which is not one"),
    fixed = TRUE)
})

test_that("a plain list still replaces the addressed field wholesale", {
  spec <- reading_spec()
  seen <- sweep_spec(spec, "fixed$coefficients", list(list(cond = 0.02)), coefficient_row)
  expect_identical(names(seen), c("coefficients", "cond"))
})

test_that("a sweep of a one-coefficient specification is unchanged", {
  spec <- build_spec(list(name = "one", seed = 5, design_kind = "between", n_subject = 40,
    factor_name = "group", lev1 = "a", lev2 = "b", intercept = 0, effect = 0.5,
    family = "gaussian", resp_name = "score", sigma = 1))
  conds <- design_conditions(effect = c(0.3, 0.6))
  merged <- sweep_spec(spec, "fixed$coefficients", conds, power_design, n_sims = 40)
  replaced <- sweep_spec(spec, "fixed$coefficients", unclass(conds), power_design, n_sims = 40)
  expect_identical(merged, replaced)
})

test_that("subsetting design_conditions() output keeps it merging, and it prints as a list", {
  spec <- reading_spec()
  conds <- design_conditions(cond = c(0.02, 0.04))
  expect_s3_class(conds, "pilotr_conditions")
  effects_only <- conds[-1]
  expect_s3_class(effects_only, "pilotr_conditions")
  expect_identical(unclass(effects_only), unclass(conds)[-1])
  seen <- sweep_spec(spec, "fixed$coefficients", effects_only, coefficient_row)
  expect_identical(names(seen), c("coefficients", names(spec$fixed$coefficients)))
  expect_identical(seen$cond, c(0.02, 0.04))
  expect_identical(utils::capture.output(print(conds)),
                   utils::capture.output(print(unclass(conds))))
})

test_that("joining design_conditions() outputs with c() keeps them merging", {
  spec <- reading_spec()
  coefs <- spec$fixed$coefficients
  first <- design_conditions(cond = 0.02)
  second <- design_conditions(age = 0.1, .null = FALSE)
  joined <- c(first, second)
  expect_s3_class(joined, "pilotr_conditions")
  expect_identical(unclass(joined), c(unclass(first), unclass(second)))
  seen <- sweep_spec(spec, "fixed$coefficients", joined, coefficient_row)
  expect_identical(names(seen), c("coefficients", names(coefs)))
  # Each condition moves only the effect it names, and the age condition leaves cond alone.
  expect_identical(seen$cond, c(0, 0.02, coefs$cond))
  expect_identical(seen$age, c(coefs$age, coefs$age, 0.1))
  expect_identical(seen$SyntaxPC, rep(coefs$SyntaxPC, 3))
  # Collapsing to an atomic vector leaves nothing to merge, so no class is claimed for it.
  expect_false(inherits(c(first, recursive = TRUE), "pilotr_conditions"))
})

test_that("sweeping one coefficient by its path keeps the others and names the column", {
  spec <- reading_spec()
  report <- function(s) data.frame(names = paste(names(s$fixed$coefficients), collapse = " "),
                                   value = s$fixed$coefficients$cond)
  seen <- sweep_spec(spec, "fixed$coefficients$cond", c(0, 0.02, 0.04), report)
  expect_identical(names(seen), c("cond", "names", "value"))
  expect_identical(seen$cond, c(0, 0.02, 0.04))
  expect_identical(seen$value, seen$cond)
  expect_identical(unique(seen$names), paste(names(spec$fixed$coefficients), collapse = " "))
})

test_that("a sweep of one coefficient solves for a minimum detectable effect", {
  spec <- build_spec(list(name = "mde", seed = 3, design_kind = "between", n_subject = 100,
    factor_name = "group", lev1 = "a", lev2 = "b", intercept = 0, effect = 0.5,
    family = "gaussian", resp_name = "score", sigma = 1))
  curve <- sweep_spec(spec, "fixed$coefficients$effect", c(0.2, 0.4, 0.6, 0.8), power_design,
                      n_sims = 200)
  s <- solve_curve(curve, target = 0.8, transform = "identity")
  expect_identical(s$x, "effect")
  # 50 subjects a group detect d = 0.57 with a power of 0.8.
  expect_lt(abs(s$value - stats::power.t.test(n = 50, sd = 1, power = 0.8)$delta), 0.05)
})
