# How units are allocated to levels and clusters, and the one allocation a specification may not
# ask for. Each between factor is assigned to blocks of its unit on its own, so two between
# factors over one unit fell into the same or overlapping blocks. A 2 x 2 between-subjects design
# over 40 subjects gave cells of 20, 0, 0 and 20, and its effects could not be estimated apart.
# The Python twin's tests/test_allocation.py runs the same cases and expects the same text, and
# tools/parity/validate_cross.py compares the two validators on them.

between_2x2 <- function(n = 40) list(
  name = "between_2x2", seed = 1,
  units = list(subject = list(n = n)),
  factors = list(
    list(name = "A", levels = c("a1", "a2"), contrasts = list(a = c(-0.5, 0.5)),
         between = "subject"),
    list(name = "B", levels = c("b1", "b2"), contrasts = list(b = c(-0.5, 0.5)),
         between = "subject")),
  fixed = list(intercept = 0, coefficients = list(a = 0.5, b = 0.3, "a:b" = 0.2)),
  random = list(),
  response = list(family = "gaussian", name = "y", sigma = 1))

# The refusal for factors between one unit, given the quoted names, "both" or "all", and the unit.
aliased <- function(names, quantifier, unit) paste0(
  "the factors ", names, " are ", quantifier, " between '", unit, "'. pilotr assigns the levels ",
  "of each between factor to blocks of ", unit, "s on its own, so the blocks of these factors ",
  "coincide or overlap, which leaves some combinations of their levels without ", unit, "s and ",
  "confounds their effects. Encode the design as one between factor whose levels are the cells, ",
  "give it the contrast columns of these factors and key each interaction as 'a:b', as ",
  "spec_from_model() in R does for a fitted pilot. Specification version 0.4 will allocate ",
  "several between factors jointly.")

refused <- function(spec, msg) expect_error(validate_spec(spec), msg, fixed = TRUE)

# One row per subject, for counting what each subject was assigned.
first_rows <- function(d) d[!duplicated(d$subject), , drop = FALSE]

# A worked encoding from the specification, which lives in the repository beside the package.
parity_case <- function(name) {
  path <- testthat::test_path("..", "..", "..", "..", "tools", "parity", "cases",
                              paste0(name, ".json"))
  testthat::skip_if_not(file.exists(path), "the parity cases are in the repository, not the package")
  load_spec(path)
}

test_that("two factors between one unit are refused, in the words the Python twin uses", {
  s <- between_2x2()
  # What the refusal prevents: every subject in A's first block is in B's first block too.
  d <- simulate_design(s, validate = FALSE)
  expect_identical(as.vector(table(d$A, d$B)), c(20L, 0L, 0L, 20L))

  refused(s, paste0(
    "the factors 'A' and 'B' are both between 'subject'. pilotr assigns the levels of each ",
    "between factor to blocks of subjects on its own, so the blocks of these factors coincide or ",
    "overlap, which leaves some combinations of their levels without subjects and confounds their ",
    "effects. Encode the design as one between factor whose levels are the cells, give it the ",
    "contrast columns of these factors and key each interaction as 'a:b', as spec_from_model() ",
    "in R does for a fitted pilot. Specification version 0.4 will allocate several between ",
    "factors jointly."))

  s2 <- load_spec(pilotr_example("between_2group_gaussian"))
  s2$factors[[2]] <- s$factors[[2]]
  refused(s2, aliased("'group' and 'B'", "both", "subject"))
})

test_that("three factors between one unit, and two between items, are refused alike", {
  s <- between_2x2(12)
  s$factors[[3]] <- list(name = "C", levels = c("c1", "c2", "c3"),
                         contrasts = list(c1 = c(-1, 1, 0), c2 = c(-1, 0, 1)), between = "subject")
  refused(s, aliased("'A', 'B' and 'C'", "all", "subject"))

  s <- load_spec(pilotr_example("crossed_mixed_rt"))
  s$factors[[2]] <- list(name = "frequency", levels = c("low", "high"),
                         contrasts = list(freq = c(-0.5, 0.5)), between = "item")
  s$factors[[3]] <- list(name = "length", levels = c("short", "long"),
                         contrasts = list(len = c(-0.5, 0.5)), between = "item")
  refused(s, aliased("'frequency' and 'length'", "both", "item"))
})

test_that("one between factor per unit stays valid", {
  # One factor between subjects and one between items: the two-list encoding relies on it.
  s <- between_2x2()
  s$units$item <- list(n = 4)
  s$factors[[2]]$between <- "item"
  expect_silent(validate_spec(s))

  # A between factor beside a within factor.
  s <- load_spec(pilotr_example("crossed_mixed_rt"))
  s$factors[[2]] <- list(name = "group", levels = c("x", "y"),
                         contrasts = list(grp = c(-0.5, 0.5)), between = "subject")
  expect_silent(validate_spec(s))
})

test_that("a between factor's blocks follow the rule the specification states", {
  # Unit u of N takes level floor((u - 1) * L / N), so 3 levels over 10 subjects give blocks of
  # 4, 3 and 3. Equal blocks need N to be a multiple of L.
  s <- load_spec(pilotr_example("between_2group_gaussian"))
  s$units$subject$n <- 10
  s$factors[[1]] <- list(name = "group", levels = c("x", "y", "z"),
                         contrasts = list(g1 = c(-1, 1, 0), g2 = c(-1, 0, 1)), between = "subject")
  s$fixed$coefficients <- list(g1 = 1)
  expect_identical(simulate_design(s)$group, rep(c("x", "y", "z"), c(4, 3, 3)))
})

test_that("a between factor and a grouping factor over one unit nest, as the specification says", {
  s <- load_spec(pilotr_example("nested_clusters"))
  s$factors <- list(list(name = "grp", levels = c("control", "treatment"),
                         contrasts = list(g = c(-0.5, 0.5)), between = "subject"))
  s$fixed$coefficients <- list(g = 0.3)
  levels_per_site <- function(s) {
    d <- first_rows(simulate_design(s))
    as.vector(tapply(d$grp, d$site, function(x) length(unique(x))))
  }
  # 12 sites, a multiple of the 2 levels: each site lies wholly in one condition.
  expect_identical(levels_per_site(s), rep(1L, 12))
  # 5 sites: the middle site straddles the two conditions.
  s$random$site$n <- 5
  expect_identical(levels_per_site(s), c(1L, 1L, 2L, 1L, 1L))
})

test_that("a 2 x 2 between design written as one factor of cells keeps every cell", {
  s <- parity_case("between_cells_2x2")
  expect_identical(as.vector(table(first_rows(simulate_design(s))$cell)), rep(10L, 4))
  s$units$subject$n <- 80
  d <- simulate_design(s)
  expect_identical(as.vector(table(first_rows(d)$cell)), rep(20L, 4))
  expect_false(anyNA(stats::coef(stats::lm(.y ~ a + b + a_b, data = model_data(s, d)))))
})

test_that("a two-list counterbalanced design shows each subject each item once", {
  s <- parity_case("two_list_counterbalanced")
  d <- simulate_design(s)
  md <- model_data(s, d)
  expect_identical(nrow(d), 24L * 18L)
  expect_identical(anyDuplicated(d[c("subject", "item")]), 0L)
  # The condition is the product column l_g. Every item appears in both conditions across
  # subjects and every subject sees both. Each list meets each condition 108 times.
  conditions <- function(by) as.vector(tapply(md$l_g, by, function(x) length(unique(x))))
  expect_identical(conditions(md$item), rep(2L, 18))
  expect_identical(conditions(md$subject), rep(2L, 24))
  expect_identical(as.vector(table(d$list, md$l_g)), rep(108L, 4))
  expect_identical(deparse(model_formula(s)),
                   ".y ~ l_g + (1 + l_g | subject) + (1 + l_g | item)")
})

test_that("randomisation within clusters puts both arms in every site", {
  s <- parity_case("within_cluster_randomised")
  d <- simulate_design(s)
  expect_identical(nrow(d), 120L)
  expect_identical(as.vector(table(d$site, model_data(s, d)$trt)), rep(5L, 24))
  # The level labels name the site each level falls in.
  expect_identical(as.integer(substr(d$arm, 2, 3)), as.integer(d$site))
})
