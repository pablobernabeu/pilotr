# A rounded response goes through R's own round(). The specification adopts its algorithm
# (src/nmath/fround.c, R 4.0.0 and later) as the rounding rule, and the Python twin ports it as
# pilotr.core.fround(). These tests pin values on which CPython's built-in round(), working on the
# exact binary value, takes the other candidate. A change to R's algorithm would therefore show
# here before it reached the cross-language contract. The inputs are written in hexadecimal, which
# every build of R parses exactly.

# near_tie() returns the design of tools/parity/cases/round_near_tie.json with the given round. Its
# sigma is far below half a unit in the last place of 1.05, so every unrounded response is the
# double nearest 1.05.
near_tie <- function(round) list(
  spec_version = "0.3", name = "round_near_tie", seed = 1,
  units = list(subject = list(n = 10)),
  fixed = list(intercept = 1.05, coefficients = setNames(list(), character(0))),
  response = list(family = "gaussian", name = "y", sigma = 1e-17, round = round))

test_that("round() takes the decimal candidate nearer in double precision, a tie going to even", {
  # In double precision, 1.1 - 1.05 and 1.05 - 1 are equal, and 10 is even. CPython's round()
  # gives 1.1 for 1.05, 4.5 for 4.55 and 394.25 for 394.255.
  expect_identical(round(0x1.0cccccccccccdp+0, 1), 1)                     # 1.05
  expect_identical(round(-0x1.0cccccccccccdp+0, 1), -1)                   # -1.05
  expect_identical(round(0x1.2333333333333p+2, 1), 0x1.2666666666666p+2)  # 4.55 to 4.6
  expect_identical(round(0x1.8a4147ae147aep+8, 2), 0x1.8a428f5c28f5cp+8)  # 394.255 to 394.26
})

test_that("a response at a decimal tie rounds by that rule", {
  expect_identical(simulate_design(near_tie(1))$y, rep(1, 10))
})

test_that("a round that is not a whole number counts as the nearest whole number of digits", {
  # Validation refuses such a value. Without validation, round() takes floor(digits + 0.5) digits,
  # and the Python twin does the same.
  s <- near_tie(2.5)
  s$fixed$intercept <- 1.23456
  expect_error(validate_spec(s), "'response.round' must be a whole number of at least 0",
               fixed = TRUE)
  three <- s
  three$response$round <- 3
  expect_identical(simulate_design(s, validate = FALSE)$y, simulate_design(three)$y)
})
