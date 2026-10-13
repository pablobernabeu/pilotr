# The specification version this release reads, and the version a specification needs. The Python
# twin's tests/test_spec_version.py runs the same cases, apart from build_spec() and
# spec_from_model(), which exist in R alone.

two_group <- function() load_spec(pilotr_example("between_2group_gaussian"))
factor_b <- function(unit = "subject") list(name = "B", levels = c("b1", "b2"),
                                            contrasts = list(b = c(-0.5, 0.5)), between = unit)

test_that("this release reads specification version 0.4 and refuses a later one", {
  expect_identical(.SPEC_VERSION, "0.4")
  s <- two_group()
  s$spec_version <- "0.4"
  expect_silent(validate_spec(s))
  s$spec_version <- "0.5"
  expect_error(validate_spec(s), paste0(
    "this specification declares spec_version 0.5, which is newer than the 0.4 this version of ",
    "pilotr understands; please upgrade pilotr"), fixed = TRUE)
})

test_that("a specification needs 0.4 for several between factors on one unit", {
  # The other 0.4 feature, a counterbalanced factor, is tested in test-counterbalance.R.
  s <- two_group()
  expect_identical(.spec_0_4_features(s), character(0))
  expect_identical(.required_spec_version(s), "0.3")
  # A factor between subjects beside one between items needs nothing new.
  s$units$item <- list(n = 4)
  s$factors[[2]] <- factor_b("item")
  expect_identical(.required_spec_version(s), "0.3")
  s$factors[[2]]$between <- "subject"
  expect_identical(.spec_0_4_features(s), "several between factors on one unit")
  expect_identical(.required_spec_version(s), "0.4")
})

test_that("the version a specification needs is never below 0.3", {
  # A specification without a 0.3 feature still gets 0.3, the version build_spec() has written
  # since 0.3.0. Releases before 0.3 ignore the field, so a lower one would gain no reader.
  s <- two_group()
  s$spec_version <- NULL
  expect_identical(.required_spec_version(s), "0.3")
  # Malformed parts are left to the validator, as for the 0.3 features.
  expect_identical(.required_spec_version(list(factors = "none")), "0.3")
  expect_identical(.required_spec_version(list(factors = list(5, list(name = "A")))), "0.3")
})

test_that("build_spec() declares 0.3, so releases that read 0.3 read what the apps export", {
  base <- list(name = "demo", seed = 1, factor_name = "group", lev1 = "a", lev2 = "b",
               n_subject = 20, n_item = 8, intercept = 0, effect = 0.5, family = "gaussian",
               resp_name = "", sigma = 1, subj_int_sd = 0.5, subj_slope_sd = 0.2, subj_corr = 0,
               item_int_sd = 0.3, item_slope_sd = 0.1, item_corr = 0)
  for (kind in c("between", "within")) for (items in c(FALSE, TRUE)) {
    spec <- build_spec(utils::modifyList(base, list(design_kind = kind, include_items = items)))
    expect_identical(names(spec)[1], "spec_version")
    expect_identical(spec$spec_version, "0.3")
    expect_silent(validate_spec(spec))
  }
  expect_match(spec_json(spec), "^\\{\\s*\"spec_version\": \"0.3\",")
})
