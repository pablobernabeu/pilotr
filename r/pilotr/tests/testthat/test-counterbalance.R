# Counterbalanced lists, new in specification version 0.4. A within factor that sets
# `counterbalance` shows each subject each item at one level, rotated across subjects and items as
# in a Latin square: subject s meets item t at combination ((s - 1) + (t - 1)) mod C of the
# counterbalanced factors' levels, read with the first-listed factor varying slowest. Until 0.4,
# every within factor was crossed with every subject-item pair, so planning a counterbalanced study
# doubled its observations per subject with two conditions. The Python twin's
# tests/test_counterbalance.py runs the same cases and expects the same text, and
# tools/parity/validate_cross.py compares the two validators on them.

# A Gaussian version of the vignettes' priming design, counterbalanced: 24 subjects and 18 items,
# the dimensions of the two-list encoding in tools/parity/cases.
cb_spec <- function(n_subject = 24, n_item = 18) list(
  spec_version = "0.4", name = "counterbalanced", seed = 2417,
  units = list(subject = list(n = n_subject), item = list(n = n_item)),
  factors = list(list(name = "cond", levels = c("related", "unrelated"),
                      contrasts = list(cond = c(-0.5, 0.5)),
                      vary_within = c("subject", "item"), counterbalance = TRUE)),
  fixed = list(intercept = 6, coefficients = list(cond = 0.06)),
  random = list(subject = list(intercept_sd = 0.12, slopes = list(cond = 0.04)),
                item = list(intercept_sd = 0.08, slopes = list(cond = 0.02))),
  response = list(family = "gaussian", name = "logRT", sigma = 0.3))

# The same design with three conditions, which the two-list encoding cannot write.
cb3_spec <- function(n_subject = 18, n_item = 12) {
  s <- cb_spec(n_subject, n_item)
  s$factors[[1]] <- list(name = "cond", levels = c("a", "b", "c"),
                         contrasts = list(c1 = c(-0.5, 0.5, 0), c2 = c(-0.5, 0, 0.5)),
                         vary_within = c("subject", "item"), counterbalance = TRUE)
  s$fixed$coefficients <- list(c1 = 0.05, c2 = 0.08)
  s$random$subject$slopes <- list(c1 = 0.04, c2 = 0.03)
  s$random$item$slopes <- list(c1 = 0.02, c2 = 0.02)
  s
}

# The level a row takes by the rule: combination ((s - 1) + (t - 1)) mod C of C levels.
rotated <- function(d, levels) levels[(d$subject - 1 + d$item - 1) %% length(levels) + 1]

refused <- function(spec, msg) expect_error(validate_spec(spec), msg, fixed = TRUE)

parity_case <- function(name) {
  path <- testthat::test_path("..", "..", "..", "..", "tools", "parity", "cases",
                              paste0(name, ".json"))
  testthat::skip_if_not(file.exists(path),
                        "the parity cases are in the repository, not the package")
  load_spec(path)
}

# The warning for unequal lists, given its lead and the counts that are not multiples.
unequal <- function(lead, subjects = NULL, items = NULL, lists) paste0(
  lead, ", and the lists are unequal.",
  if (!is.null(subjects)) sprintf(paste0(
    " The %d subjects are not a multiple of %d, so the lists differ in size by one subject, and ",
    "an item is seen in some conditions by one subject more than in others."), subjects, lists),
  if (!is.null(items)) sprintf(paste0(
    " The %d items are not a multiple of %d, so each list shows some conditions on one item more ",
    "than others."), items, lists))

test_that("a counterbalanced factor gives each subject-item pair one row", {
  s <- cb_spec()
  expect_silent(validate_spec(s))
  d <- simulate_design(s)
  expect_identical(nrow(d), 24L * 18L)
  expect_identical(anyDuplicated(d[c("subject", "item")]), 0L)
  # The canonical order: subjects slowest, then items.
  expect_identical(d$subject, rep(1:24, each = 18))
  expect_identical(d$item, rep(1:18, times = 24))
  expect_identical(names(d), c("subject", "item", "cond", "logRT"))
})

test_that("subject s meets item t at level ((s - 1) + (t - 1)) mod C", {
  d <- simulate_design(cb_spec())
  expect_identical(d$cond, rotated(d, c("related", "unrelated")))
  d3 <- simulate_design(cb3_spec())
  expect_identical(d3$cond, rotated(d3, c("a", "b", "c")))
  # Subject 1 follows the first list, item 1 at the first level, and subject 2 the second list.
  expect_identical(d3$cond[d3$subject == 1][1:4], c("a", "b", "c", "a"))
  expect_identical(d3$cond[d3$subject == 2][1:4], c("b", "c", "a", "b"))
})

test_that("every item appears in every condition, and every subject sees each equally often", {
  # 18 subjects and 12 items, both multiples of the three lists.
  d <- simulate_design(cb3_spec())
  expect_identical(as.vector(table(d$item, d$cond)), rep(6L, 12 * 3))
  expect_identical(as.vector(table(d$subject, d$cond)), rep(4L, 18 * 3))
  # With 14 items, a subject sees some conditions on one item more than others.
  s <- cb3_spec(18, 14)
  d <- suppressWarnings(simulate_design(s))
  per_subject <- table(d$subject, d$cond)
  expect_identical(sort(unique(as.vector(per_subject))), c(4L, 5L))
  expect_identical(as.vector(rowSums(per_subject)), rep(14, 18))
})

test_that("several counterbalanced factors rotate jointly, the first-listed slowest", {
  s <- cb_spec(8, 8)
  s$factors[[2]] <- list(name = "prime", levels = c("p1", "p2"),
                         contrasts = list(pr = c(-0.5, 0.5)),
                         vary_within = c("subject", "item"), counterbalance = TRUE)
  s$fixed$coefficients$pr <- 0.03
  expect_silent(validate_spec(s))
  d <- simulate_design(s)
  expect_identical(nrow(d), 64L)
  cell <- (d$subject - 1 + d$item - 1) %% 4
  expect_identical(d$cond, c("related", "unrelated")[cell %/% 2 + 1])
  expect_identical(d$prime, c("p1", "p2")[cell %% 2 + 1])
})

test_that("beside a counterbalanced factor, only the crossed within factors expand rows", {
  s <- cb_spec(4, 6)
  s$factors[[2]] <- list(name = "task", levels = c("read", "listen"),
                         contrasts = list(tk = c(-0.5, 0.5)), vary_within = c("subject", "item"))
  s$fixed$coefficients$tk <- 0.1
  d <- simulate_design(s)
  expect_identical(nrow(d), 4L * 6L * 2L)
  # Each subject-item pair has one row per level of the crossed factor, in level order, and keeps
  # the one condition its rotation gives it.
  expect_identical(d$subject, rep(1:4, each = 12))
  expect_identical(d$item, rep(rep(1:6, each = 2), times = 4))
  expect_identical(d$task, rep(c("read", "listen"), times = 24))
  expect_identical(d$cond, rotated(d, c("related", "unrelated")))
})

test_that("counterbalancing draws nothing, so it gives the responses of an encoding by unit", {
  # The rotation of two conditions, written with a factor between subjects and one between items
  # whose contrasts multiply to the condition, each subject and each item a level of its own. The
  # two specifications have the same rows in the same order, the same draws and the same linear
  # predictor term by term, so their responses agree to the last bit.
  native <- cb_spec()
  encoded <- native
  encoded$spec_version <- "0.3"
  encoded$factors <- list(
    list(name = "subject_list", levels = sprintf("s%02d", 1:24),
         contrasts = list(l = rep(c(1, -1), 12)), between = "subject"),
    list(name = "item_list", levels = sprintf("t%02d", 1:18),
         contrasts = list(g = rep(c(-0.5, 0.5), 9)), between = "item"))
  encoded$fixed$coefficients <- list(`l:g` = 0.06)
  encoded$random$subject$slopes <- list(`l:g` = 0.04)
  encoded$random$item$slopes <- list(`l:g` = 0.02)
  expect_identical(simulate_design(native)$logRT, simulate_design(encoded)$logRT)
})

test_that("two counterbalanced conditions have the cell counts of the two-list encoding", {
  native <- simulate_design(parity_case("counterbalanced_two_conditions"))
  two_list <- parity_case("two_list_counterbalanced")
  encoded <- model_data(two_list, simulate_design(two_list))
  expect_identical(nrow(native), nrow(encoded))
  counts <- function(unit, condition) as.vector(table(unit, condition))
  # Nine items in each condition for every subject, and 12 subjects in each for every item.
  expect_identical(counts(native$subject, native$cond), counts(encoded$subject, encoded$l_g))
  expect_identical(counts(native$item, native$cond), counts(encoded$item, encoded$l_g))
  # Subject s follows list (s - 1) mod 2, and each list meets each condition 108 times.
  expect_identical(counts((native$subject - 1) %% 2, native$cond),
                   counts(encoded$list, encoded$l_g))
})

test_that("a counterbalanced design is analysed with the model of the crossed design", {
  # Each subject meets every condition across items and each item every condition across
  # subjects, so the by-subject and by-item slopes on the condition stay in the model, and
  # model_formula() and model_data() read the specification as they read the crossed one.
  cb <- cb_spec()
  crossed <- cb
  crossed$factors[[1]]$counterbalance <- NULL
  expect_identical(model_formula(cb), model_formula(crossed))
  expect_identical(deparse(model_formula(cb)),
                   ".y ~ cond + (1 + cond || subject) + (1 + cond || item)")
  expect_identical(names(model_data(cb, simulate_design(cb))),
                   names(model_data(crossed, simulate_design(crossed))))
})

test_that("the parity cases hold the designs they describe", {
  d <- simulate_design(parity_case("counterbalanced_two_conditions"))
  expect_identical(nrow(d), 432L)
  expect_identical(anyDuplicated(d[c("subject", "item")]), 0L)
  expect_identical(d$cond, rotated(d, c("related", "unrelated")))

  d <- simulate_design(parity_case("counterbalanced_three_conditions"))
  expect_identical(nrow(d), 18L * 12L * 2L)
  expect_identical(d$cond, rotated(d, c("a", "b", "c")))
  expect_identical(d$task, rep(c("read", "listen"), times = 18 * 12))
  expect_identical(as.vector(table(d$item, d$cond)), rep(12L, 12 * 3))
})

test_that("counterbalance false and an absent field cross the factor alike", {
  s <- cb_spec()
  s$factors[[1]]$counterbalance <- FALSE
  expect_silent(validate_spec(s))
  crossed <- s
  crossed$factors[[1]]$counterbalance <- NULL
  expect_identical(simulate_design(s), simulate_design(crossed))
  expect_identical(nrow(simulate_design(s)), 24L * 18L * 2L)
})

test_that("counterbalance needs specification version 0.4, whatever its value", {
  msg <- paste0("this specification uses factors counterbalance, which requires spec_version ",
                "\"0.4\", but declares 0.3. An implementation that reads no version above 0.3 ",
                "refuses such a field or ignores it, and ignoring it generates different data.")
  s <- cb_spec()
  s$spec_version <- "0.3"
  refused(s, msg)
  expect_error(validate_spec(s, strict = FALSE), msg, fixed = TRUE)
  # Releases that read 0.3 do not know the field, whatever its value.
  s$factors[[1]]$counterbalance <- FALSE
  refused(s, msg)
  s$factors[[1]]$counterbalance <- TRUE
  s$spec_version <- NULL
  refused(s, sub("declares 0.3", "declares 0.2", msg, fixed = TRUE))
})

test_that("counterbalance applies to a within factor of a fully crossed design with items", {
  s <- cb_spec()
  s$factors[[2]] <- list(name = "group", levels = c("x", "y"), contrasts = list(grp = c(-0.5, 0.5)),
                         between = "subject", counterbalance = TRUE)
  refused(s, "factors[2].counterbalance applies only to a within factor, one that sets vary_within")

  s <- cb_spec()
  s$units$item <- NULL
  s$random$item <- NULL
  s$factors[[1]]$vary_within <- "subject"
  refused(s, paste0("factors[1].counterbalance rotates the levels over items, but the design has ",
                    "no item unit"))

  s <- cb_spec()
  s$units$item$per_subject <- 6
  refused(s, paste0("factors[1].counterbalance needs every subject to see every item, so it ",
                    "cannot be combined with units.item.per_subject"))

  s <- cb_spec()
  s$factors[[1]]$counterbalance <- "yes"
  refused(s, "factors[1].counterbalance must be true or false")
  s$factors[[1]]$counterbalance <- NA
  refused(s, "factors[1].counterbalance must be true or false")
})

test_that("a design whose subjects or items are not a multiple of the lists warns", {
  lead <- "the counterbalanced factor 'cond' rotates its 3 levels over 3 lists"
  s <- cb3_spec(18, 12)
  expect_silent(validate_spec(s))
  s$units$subject$n <- 20
  expect_warning(validate_spec(s), unequal(lead, subjects = 20, lists = 3), fixed = TRUE)
  s$units$subject$n <- 18
  s$units$item$n <- 16
  expect_warning(validate_spec(s), unequal(lead, items = 16, lists = 3), fixed = TRUE)
  s$units$subject$n <- 20
  expect_warning(validate_spec(s, strict = FALSE), unequal(lead, 20, 16, 3), fixed = TRUE)
  # The data are still the rotation, with lists one subject or one item apart.
  d <- suppressWarnings(simulate_design(s))
  expect_identical(d$cond, rotated(d, c("a", "b", "c")))

  # Two counterbalanced factors rotate the combinations of their levels.
  s <- cb_spec(6, 8)
  s$factors[[2]] <- list(name = "prime", levels = c("p1", "p2"),
                         contrasts = list(pr = c(-0.5, 0.5)),
                         vary_within = c("subject", "item"), counterbalance = TRUE)
  expect_warning(validate_spec(s), unequal(paste0(
    "the counterbalanced factors 'cond' and 'prime' rotate the 4 combinations of their levels ",
    "over 4 lists"), subjects = 6, lists = 4), fixed = TRUE)
})

test_that("the version a counterbalanced specification needs is 0.4", {
  s <- cb_spec()
  expect_identical(.spec_0_4_features(s), "factors counterbalance")
  expect_identical(.required_spec_version(s), "0.4")
  s$factors[[1]]$counterbalance <- FALSE
  expect_identical(.required_spec_version(s), "0.4")
  s$factors[[1]]$counterbalance <- NULL
  expect_identical(.spec_0_4_features(s), character(0))
  expect_identical(.required_spec_version(s), "0.3")
})

test_that("a counterbalance written as jsonlite writes it reads as its value", {
  # jsonlite::write_json() writes every scalar as a one-element array by default, and the Python
  # twin reads [true] as true, as R does. A null means the field is absent.
  boxed <- sub("\"counterbalance\": true", "\"counterbalance\": [true]", spec_json(cb_spec()),
               fixed = TRUE)
  expect_true(grepl("[true]", boxed, fixed = TRUE))
  expect_identical(simulate_design(load_spec(boxed))$logRT, simulate_design(cb_spec())$logRT)
  absent <- sub("\"counterbalance\": true", "\"counterbalance\": null", spec_json(cb_spec()),
                fixed = TRUE)
  expect_identical(nrow(simulate_design(load_spec(absent))), 24L * 18L * 2L)
})

test_that("build_spec() offers counterbalanced lists for a within design with items", {
  p <- list(name = "cb", seed = 1, design_kind = "within", include_items = TRUE,
            n_subject = 12, n_item = 8, factor_name = "cond", lev1 = "a", lev2 = "b",
            intercept = 0, effect = 0.5, family = "gaussian", resp_name = "", sigma = 1,
            subj_int_sd = 0.5, subj_slope_sd = 0.2, subj_corr = 0, item_int_sd = 0.3,
            item_slope_sd = 0.1, item_corr = 0, counterbalance = TRUE)
  spec <- build_spec(p)
  expect_identical(spec$spec_version, "0.4")
  expect_true(spec$factors[[1]]$counterbalance)
  expect_silent(validate_spec(spec))
  d <- simulate_design(spec)
  expect_identical(nrow(d), 12L * 8L)
  expect_identical(d$cond, rotated(d, c("a", "b")))
  expect_match(spec_json(spec), "\"counterbalance\": true", fixed = TRUE)
  expect_identical(simulate_design(load_spec(spec_json(spec))), d)
  # Where the option cannot apply it is ignored, as when an app hides a ticked box after the
  # design changes. The specification then needs nothing above 0.3.
  for (change in list(list(design_kind = "between"), list(include_items = FALSE))) {
    spec <- build_spec(utils::modifyList(p, change))
    expect_null(spec$factors[[1]]$counterbalance)
    expect_identical(spec$spec_version, "0.3")
    expect_silent(validate_spec(spec))
  }
})
