# Names, columns and structure. Each specification below used to validate and then move an
# effect, rescale a variance or overwrite a column without a word. The Python twin's
# tests/test_validate_names.py runs the same cases and expects the same text, and
# tools/parity/validate_cross.py compares the two validators on them.

names_between <- function() load_spec(pilotr_example("between_2group_gaussian"))
names_crossed <- function() load_spec(pilotr_example("crossed_mixed_rt"))

COLUMN_OF_ITS_OWN <- paste0("; every unit, grouping factor, factor, predictor and the response ",
                            "needs a column of its own")
OWN_FACTOR_ONLY <- "; a contrast column may share its own factor's name but no other column's"

refused <- function(spec, msg) expect_error(validate_spec(spec), msg, fixed = TRUE)

test_that("two columns with one name are refused", {
  s <- names_between()
  s$factors[[1]]$name <- "subject"
  refused(s, paste0("the name 'subject' is used for more than one column (the subject unit and ",
                    "factors[1])", COLUMN_OF_ITS_OWN))

  s <- names_between()
  s$response$name <- "group"
  refused(s, paste0("the name 'group' is used for more than one column (factors[1] and the ",
                    "response)", COLUMN_OF_ITS_OWN))

  s <- load_spec(pilotr_example("nested_clusters"))
  s$predictors <- list(list(name = "site", varies_by = "subject"))
  s$response$name <- "site"
  refused(s, paste0("the name 'site' is used for more than one column (random.site, ",
                    "predictors[1] and the response)", COLUMN_OF_ITS_OWN))
})

test_that("a contrast column belongs to one factor and shares no other column's name", {
  s <- names_between()
  s$factors[[1]]$contrasts <- list(x = c(-0.5, 0.5))
  s$factors[[2]] <- list(name = "block", levels = c("p", "q"),
                         contrasts = list(x = c(-0.5, 0.5)), between = "subject")
  s$fixed$coefficients <- list(x = 5)
  refused(s, paste0("the contrast column 'x' is defined by more than one factor (factors[1] and ",
                    "factors[2]); each contrast column has to belong to a single factor"))

  # A within factor reusing the between factor's column put the whole effect on the block.
  s <- names_between()
  s$factors[[1]]$contrasts <- list(effect = c(-0.5, 0.5))
  s$factors[[2]] <- list(name = "block", levels = c("p", "q"),
                         contrasts = list(effect = c(-0.5, 0.5)), vary_within = "subject")
  s$fixed$coefficients <- list(effect = 5)
  refused(s, paste0("the contrast column 'effect' is defined by more than one factor ",
                    "(factors[1] and factors[2]); each contrast column has to belong to a ",
                    "single factor"))

  # A predictor named like the contrast replaced the group difference with its own draws.
  s <- names_between()
  s$factors[[1]]$contrasts <- list(effect = c(-0.5, 0.5))
  s$predictors <- list(list(name = "effect", varies_by = "subject"))
  s$fixed$coefficients <- list(effect = 5)
  refused(s, paste0("the contrast column 'effect' of factors[1] is also the name of ",
                    "predictors[1]", OWN_FACTOR_ONLY))

  # A contrast named 'item' replaced the item identifiers in model_data(), so the formula's
  # (1 | item) grouped by condition.
  s <- names_crossed()
  s$factors[[1]]$contrasts$item <- c(-0.5, 0.5)
  refused(s, paste0("the contrast column 'item' of factors[1] is also the name of the item unit",
                    OWN_FACTOR_ONLY))

  s <- names_between()
  s$factors[[2]] <- list(name = "block", levels = c("p", "q"),
                         contrasts = list(group = c(-0.5, 0.5)), between = "subject")
  refused(s, paste0("the contrast column 'group' of factors[2] is also the name of factors[1]",
                    OWN_FACTOR_ONLY))

  # A contrast named after its own factor stays allowed: two-level specifications use it.
  s <- names_between()
  s$factors[[1]]$contrasts <- list(group = c(-0.5, 0.5))
  s$fixed$coefficients <- list(group = 5)
  expect_silent(validate_spec(s))
})

test_that("an interaction's analysis column may not take another column's name", {
  s <- names_crossed()
  s$predictors <- list(list(name = "freq", varies_by = "item"),
                       list(name = "cond_freq", varies_by = "subject"))
  s$fixed$coefficients[["cond:freq"]] <- 0.01
  refused(s, paste0("the interaction 'cond:freq' becomes the analysis column 'cond_freq', which ",
                    "is already the name of predictors[2]"))

  s <- names_crossed()
  s$factors[[1]]$contrasts$cond_cond2 <- c(1, -1)
  s$factors[[1]]$contrasts$cond2 <- c(-1, 1)
  s$fixed$coefficients[["cond:cond2"]] <- 0.01
  refused(s, paste0("the interaction 'cond:cond2' becomes the analysis column 'cond_cond2', ",
                    "which is already the name of a contrast column of factors[1]"))
})

test_that("a level listed twice is refused", {
  s <- names_between()
  s$factors[[1]]$levels <- c("x", "x")
  refused(s, "factors[1].levels repeats 'x'; each level needs its own label")
})

test_that("a factor both between and within a unit is refused", {
  s <- names_crossed()
  s$factors[[1]]$between <- "subject"
  refused(s, "factors[1] sets both 'vary_within' and 'between'; a factor has to set exactly one of them")
})

test_that("a correlation pairs two different terms, and each pair once", {
  s <- names_crossed()
  s$random$subject$correlations <- list(`cond,cond` = 0.25)
  refused(s, paste0("random.subject.correlations key 'cond,cond' pairs 'cond' with itself; a ",
                    "term's correlation with itself is always 1"))

  s <- names_crossed()
  s$random$subject$correlations <- list(`intercept,cond` = 0.9, `cond~intercept` = -0.9)
  refused(s, paste0("random.subject.correlations keys 'intercept,cond' and 'cond~intercept' ",
                    "name the same pair of terms; give each pair once"))
})

test_that("random.item is refused in a design without items", {
  s <- names_between()
  s$random <- list(item = list(intercept_sd = 0.5))
  refused(s, paste0("random.item describes an item unit the design does not have; add ",
                    "units.item or remove random.item"))
})

test_that("a blank name is refused where it used to fail inside the simulator", {
  # Ported from 7693ab2 on claude/meridian-packages-apps-review-g2a74i. An emptied name built a
  # spec that validated and then stopped in simulate_design() with base R's "replacement has
  # length zero", where the twin wrote a column with no name at all.
  s <- names_between()
  for (blank in c("", "  ")) {
    s$factors[[1]]$name <- blank
    refused(s, "factors[1].name must be a non-empty string")
  }
  s <- load_spec(pilotr_example("reading_time_continuous"))
  s$predictors[[1]]$name <- ""
  refused(s, "predictors[1].name must be a non-empty string")
  s <- names_between()
  s$response$name <- " "
  refused(s, "'response.name' must be a non-empty string")
  s <- load_spec(pilotr_example("nested_clusters"))
  names(s$random)[names(s$random) == "site"] <- ""
  refused(s, "a 'random' grouping factor must have a non-empty name")
})

test_that("two messages read as the twin's do", {
  s <- names_between()
  s$response <- list(family = "ordinal", name = "r", thresholds = "low")
  refused(s, "'response.thresholds' must be a number or a non-empty numeric array")
})

test_that("a within factor that omits a unit of the design warns, in either mode", {
  s <- names_crossed()
  full <- simulate_design(s)
  s$factors[[1]]$vary_within <- "subject"
  msg <- paste0("factors[1].vary_within lists 'subject' but not 'item'. pilotr crosses a within ",
                "factor with every unit of the design, so it varies within 'item' as well; list ",
                "every unit, or make it between 'item' if items carry it. From spec_version 0.4 ",
                "this is an error.")
  expect_warning(validate_spec(s), msg, fixed = TRUE)
  expect_warning(validate_spec(s, strict = FALSE), msg, fixed = TRUE)
  # The list never changed the data, which is why it is a warning and not yet an error.
  expect_identical(suppressWarnings(simulate_design(s)), full)

  # A list naming every unit, in any order, and 'subject' alone in a design without items are
  # complete.
  s$factors[[1]]$vary_within <- c("item", "subject")
  expect_silent(validate_spec(s))
  s <- names_between()
  s$factors[[1]]$between <- NULL
  s$factors[[1]]$vary_within <- "subject"
  expect_silent(validate_spec(s))
})

test_that("model_data() reads every factor's levels before writing any contrast column", {
  # The first contrast column carried the factor's own name and overwrote the labels, so the
  # second matched numbers against labels and came out NA in every row.
  s <- list(name = "three", seed = 1,
            units = list(subject = list(n = 6L), item = list(n = 3L)),
            factors = list(list(name = "cond", levels = c("a", "b", "c"),
                                contrasts = list(cond = c(-1, 1, 0), cond2 = c(-1, 0, 1)),
                                vary_within = c("subject", "item"))),
            fixed = list(intercept = 0, coefficients = list(cond = 0.5, cond2 = 0.2)),
            random = list(),
            response = list(family = "gaussian", name = "y", sigma = 1))
  md <- model_data(s, simulate_design(s))
  expect_false(anyNA(md$cond))
  expect_false(anyNA(md$cond2))
  expect_setequal(unique(md$cond2), c(-1, 0, 1))

  # The same overwrite across factors, which validate_spec() now refuses but model_data() may
  # still be handed directly.
  s2 <- s
  s2$factors <- list(list(name = "A", levels = c("a1", "a2"), contrasts = list(B = c(-0.5, 0.5)),
                          vary_within = c("subject", "item")),
                     list(name = "B", levels = c("b1", "b2"), contrasts = list(bc = c(-0.5, 0.5)),
                          between = "subject"))
  s2$fixed$coefficients <- list(B = 0.5, bc = 0.2)
  md2 <- model_data(s2, simulate_design(s2, validate = FALSE))
  expect_false(anyNA(md2$bc))
})
