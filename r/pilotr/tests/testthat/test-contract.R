# The cross-language contract as the documentation states it, and the schema against strict
# validation. The guarantee was stated as bit-identity for the Gaussian family and for every
# family with rounding set, and neither held. AS 241's tail branch passes about 15% of normal
# draws through log(), whose last bit each platform's maths library decides, and validate_spec()
# refuses round for the discrete families. The schema accepted an unknown top-level field and a
# response parameter its family does not use, which strict validation refuses, and
# ?replicate_seeds gave the wrong reason why consecutive seeds correlate. The checks on files
# outside the installed package skip when the repository is not there. Twinned with
# python/tests/test_schema.py and python/tests/test_contract.py.

repo_file <- function(...) {
  path <- testthat::test_path("..", "..", ...)
  testthat::skip_if_not(file.exists(path), paste(file.path(...), "not available"))
  path
}

# One line of text, so that a sentence wrapped anywhere is still found.
flat_text <- function(path) {
  gsub("\\s+", " ", paste(readLines(path, encoding = "UTF-8", warn = FALSE), collapse = " "))
}

schema <- function() {
  path <- repo_file("..", "..", "spec", "design.schema.json")
  jsonlite::fromJSON(path, simplifyVector = FALSE)
}

# Each family's required and refused response fields, read off the schema's if/then blocks.
family_rules <- function(s) {
  rules <- list()
  for (block in s[["properties"]][["response"]][["allOf"]]) {
    cond <- block[["if"]][["properties"]][["family"]]
    fams <- if (!is.null(cond[["const"]])) cond[["const"]] else unlist(cond[["enum"]])
    then <- block[["then"]]
    for (fam in fams) {
      r <- rules[[fam]]
      if (is.null(r)) r <- list(required = character(0), refused = character(0), round_min = 0)
      r$required <- c(r$required, unlist(then[["required"]]))
      for (key in names(then[["properties"]])) {
        sub <- then[["properties"]][[key]]
        if (isFALSE(sub)) r$refused <- c(r$refused, key)
        else if (identical(key, "round")) r$round_min <- sub[["minimum"]]
      }
      rules[[fam]] <- r
    }
  }
  rules
}

test_that("the schema refuses unknown top-level fields, as strict validation does", {
  expect_false(schema()[["additionalProperties"]])
})

test_that("the schema requires and refuses each family's parameters as validate_spec() does", {
  rules <- family_rules(schema())
  family_params <- pilotr:::.family_params
  expect_setequal(names(rules), names(family_params))
  params <- c("sigma", "shift", "beta", "phi", "thresholds")
  for (fam in names(family_params)) {
    needed <- family_params[[fam]]
    refused <- setdiff(params, needed)
    if (!fam %in% pilotr:::.rounding_families) refused <- c(refused, "round")
    expect_setequal(rules[[fam]]$required, needed)
    expect_setequal(rules[[fam]]$refused, refused)
  }
  # A proportion rounded to whole numbers is 0 or 1, outside the family's support.
  expect_equal(rules[["beta"]]$round_min, 1)
})

test_that("the schema accepts the single-value forms strict validation accepts", {
  s <- schema()
  branch_types <- function(node) vapply(node[["anyOf"]], function(b) b[["type"]], character(1))
  vary_within <- s[["properties"]][["factors"]][["items"]][["properties"]][["vary_within"]]
  expect_setequal(branch_types(vary_within), c("array", "string"))
  thresholds <- s[["properties"]][["response"]][["properties"]][["thresholds"]]
  expect_setequal(branch_types(thresholds), c("array", "number"))
})

test_that("DESCRIPTION states the guarantee with its platform condition", {
  desc <- read.dcf(system.file("DESCRIPTION", package = "pilotr"), fields = "Description")
  desc <- gsub("\\s+", " ", desc[1, 1])
  expect_false(grepl("for any family with rounding set", desc, fixed = TRUE))
  expect_true(grepl("maths library", desc, fixed = TRUE))
})

test_that("the vignettes show determinism with identical(), not all.equal()", {
  # all.equal() allows a relative difference of 1.5e-8, so it cannot show that two runs agree.
  for (v in c("cross-language.Rmd", "getting-started.Rmd")) {
    text <- flat_text(repo_file("vignettes", v))
    expect_false(grepl("all.equal(simulate_design", text, fixed = TRUE), info = v)
    expect_true(grepl("identical(simulate_design(spec), simulate_design(spec))", text,
                      fixed = TRUE), info = v)
  }
  expect_false(grepl("agree to full double precision",
                     flat_text(repo_file("vignettes", "cross-language.Rmd")), fixed = TRUE))
  # The within-language guarantee holds on one platform at a time.
  expect_false(grepl("in Python and on any machine",
                     flat_text(repo_file("vignettes", "getting-started.Rmd")), fixed = TRUE))
})

test_that("?replicate_seeds gives the affine mechanism", {
  text <- flat_text(repo_file("man", "replicate_seeds.Rd"))
  expect_false(grepl("a few steps apart", text, fixed = TRUE))
  expect_true(grepl("affine", text, fixed = TRUE))
  expect_true(grepl("10.1145/1276927.1276928", text, fixed = TRUE))
})

test_that("?load_spec and ?simulate_design say what any other value of validate does", {
  # A value other than TRUE or FALSE validates leniently in both functions. ?load_spec said that
  # it was passed to validate_spec() as its `strict` argument, and ?simulate_design said nothing.
  spec <- load_spec(pilotr_example("between_2group_gaussian"))
  spec$extra_field <- "x"
  tmp <- tempfile(fileext = ".json")
  on.exit(unlink(tmp))
  writeLines(spec_json(spec), tmp)
  expect_error(load_spec(tmp), "unknown top-level field 'extra_field'")
  for (v in list(NA, 1, "lenient"))
    expect_warning(load_spec(tmp, validate = v), "unknown top-level field 'extra_field'")
  expect_warning(simulate_design(spec, validate = NA), "unknown top-level field 'extra_field'")
  expect_silent(load_spec(tmp, validate = FALSE))
  for (topic in c("load_spec.Rd", "simulate_design.Rd")) {
    text <- flat_text(repo_file("man", topic))
    expect_false(grepl("as its \\code{strict} argument", text, fixed = TRUE), info = topic)
    expect_true(grepl("Any other value validates leniently", text, fixed = TRUE), info = topic)
  }
})

test_that("?pilotr names every response family", {
  # The overview left out the ex-Gaussian, new in 0.3.
  text <- flat_text(repo_file("man", "pilotr-package.Rd"))
  for (family in c("Gaussian", "lognormal", "shifted lognormal", "ex-Gaussian", "Bernoulli",
                   "Poisson", "ordinal", "Beta"))
    expect_true(grepl(family, text, fixed = TRUE), info = family)
  expect_length(pilotr:::.family_params, 8L)
})

test_that("R's ^ computes a square as u * u, as core.py's boost for a shape of 1/2 does", {
  # core.py's gamma_mt() computes the boost u^(1 / shape) as u * u when the exponent is exactly 2,
  # since R's ^ does, where Python's ** calls pow(). glibc's pow() returned another double for
  # 0.08% of this generator's uniforms. Should R stop squaring by multiplication, core.py and the
  # paragraph of spec/SPEC.md under Scope of the guarantee need revising.
  rng <- make_rng(1)
  u <- vapply(1:2000, function(i) rng$uniform(), numeric(1))
  expect_identical(u^(1 / 0.5), u * u)
  # The 860th uniform of seed 1, for which glibc's pow() returned the double below u * u. With
  # normals of 0, .gamma_mt() accepts its first proposal, 3/2 - 1/3, and boosts it with this
  # uniform, as test_basics.py checks for core.py.
  expect_identical(u[860], 0.6979060118654794)
  draws <- list(normal = function() 0, uniform = function() u[860])
  expect_identical(.gamma_mt(draws, 0.5), (1.5 - 1 / 3) * (u[860] * u[860]))
})

test_that("R's ^ computes a cube and a fourth power as the specification says", {
  # spec/SPEC.md, under Scope of the guarantee, says that from R 4.6.0 `^` computes a cube or a
  # fourth power by repeated multiplication, where Python's ** calls pow(), so a Beta shape of
  # exactly 1/3 or 1/4 can differ between the twins on any platform. Should R change this again,
  # that paragraph needs revising.
  skip_if(getRversion() < "4.6.0", "R before 4.6.0 calls the maths library for these powers")
  rng <- make_rng(4242)
  u <- vapply(1:1000, function(i) rng$uniform(), numeric(1))
  expect_identical(u^3, u * u * u)
  expect_identical(u^4, u * u * u * u)
})

test_that("consecutive seeds step by a constant, as ?replicate_seeds says", {
  lag1 <- function(x) stats::cor(x[-length(x)], x[-1])
  first <- vapply(1000:2999, function(s) make_rng(s)$uniform(), numeric(1))
  steps <- (diff(first) + 0.5) %% 1 - 0.5
  # The step is the same for every pair of consecutive seeds, to within 1e-6.
  expect_lt(max(abs(steps + 0.009056)), 1e-6)
  d <- abs(mean(steps))
  expect_equal(lag1(first), 1 - 6 * d * (1 - d), tolerance = 0.01)
  expect_gt(lag1(first), 0.94)
  # Seeds drawn from the generator do not step, so their first draws barely correlate.
  seeds <- replicate_seeds(90210, 2000)
  expect_lt(abs(lag1(vapply(seeds, function(s) make_rng(s)$uniform(), numeric(1)))), 0.05)
})
