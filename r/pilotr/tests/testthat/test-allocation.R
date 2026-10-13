# How units are allocated to levels and clusters. Up to specification version 0.3, each between
# factor was assigned to blocks of its unit on its own, so two between factors over one unit fell
# into the same or overlapping blocks. A 2 x 2 between-subjects design over 40 subjects gave cells
# of 20, 0, 0 and 20, and its effects could not be estimated apart, so a specification declaring
# a version below 0.4 is refused for them. Version 0.4 allocates them jointly, one block of units
# per combination of their levels. The Python twin's tests/test_allocation.py runs the same cases
# and expects the same text, and tools/parity/validate_cross.py compares the two validators on
# them.

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
  "the factors ", names, " are ", quantifier, " between '", unit, "'. Below specification ",
  "version 0.4, pilotr assigns the levels of each between factor to blocks of ", unit, "s on its ",
  "own. The blocks of these factors then coincide or overlap, which leaves some combinations of ",
  "their levels without ", unit, "s and confounds their effects. Declare spec_version \"0.4\", ",
  "which allocates several between factors jointly. Alternatively, encode the design as one ",
  "between factor whose levels are the cells, give it the contrast columns of these factors and ",
  "key each interaction as 'a:b', as spec_from_model() in R does for a fitted pilot. Earlier ",
  "releases of pilotr read that encoding too.")

# The labels of each row's levels of A and B, for reading off a joint allocation.
cells_of <- function(d, sep = " ") paste(d$A, d$B, sep = sep)

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

test_that("two factors between one unit are refused below 0.4, in the words the Python twin uses", {
  s <- between_2x2()
  refused(s, paste0(
    "the factors 'A' and 'B' are both between 'subject'. Below specification version 0.4, pilotr ",
    "assigns the levels of each between factor to blocks of subjects on its own. The blocks of ",
    "these factors then coincide or overlap, which leaves some combinations of their levels ",
    "without subjects and confounds their effects. Declare spec_version \"0.4\", which allocates ",
    "several between factors jointly. Alternatively, encode the design as one between factor ",
    "whose levels are the cells, give it the contrast columns of these factors and key each ",
    "interaction as 'a:b', as spec_from_model() in R does for a fitted pilot. Earlier releases of ",
    "pilotr read that encoding too."))
  # Declaring 0.3 changes nothing: that version's rule put these factors in the same blocks,
  # leaving cells of 20, 0, 0 and 20 subjects.
  s$spec_version <- "0.3"
  refused(s, aliased("'A' and 'B'", "both", "subject"))

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

# Version 0.4 allocates the between factors of one unit jointly. The units fall into one block per
# combination of levels, a cell, in order, with the first-listed factor varying slowest.

test_that("a 0.4 specification allocates two factors between subjects jointly", {
  s <- between_2x2()
  s$spec_version <- "0.4"
  expect_silent(validate_spec(s))
  d <- simulate_design(s)
  # Subjects 1 to 10 take the first cell, 11 to 20 the second and so on, B varying fastest.
  expect_identical(cells_of(d), rep(c("a1 b1", "a1 b2", "a2 b1", "a2 b2"), each = 10))
  # Both effects and the interaction can be estimated.
  expect_false(anyNA(stats::coef(stats::lm(.y ~ a + b + a_b, data = model_data(s, d)))))
})

test_that("a 3 x 2 design over 12 subjects puts two subjects in each cell", {
  s <- between_2x2(12)
  s$spec_version <- "0.4"
  s$factors[[1]] <- list(name = "A", levels = c("a1", "a2", "a3"),
                         contrasts = list(a1 = c(-0.5, 0.5, 0), a2 = c(-0.5, 0, 0.5)),
                         between = "subject")
  s$fixed$coefficients <- list(a1 = 0.4, a2 = 0.7, b = 0.3)
  expect_identical(cells_of(simulate_design(s)),
                   rep(c("a1 b1", "a1 b2", "a2 b1", "a2 b2", "a3 b1", "a3 b2"), each = 2))
})

test_that("the cells are equal only when the units are a multiple of the cells", {
  # Unit u of 10 takes cell floor((u - 1) * 4 / 10), so the four cells hold 3, 2, 3 and 2.
  s <- between_2x2(10)
  s$spec_version <- "0.4"
  expect_identical(as.vector(table(cells_of(simulate_design(s)))), c(3L, 2L, 3L, 2L))
})

test_that("two factors between items are allocated jointly too", {
  s <- load_spec(pilotr_example("crossed_mixed_rt"))
  s$spec_version <- "0.4"
  s$factors[[2]] <- list(name = "frequency", levels = c("low", "high"),
                         contrasts = list(freq = c(-0.5, 0.5)), between = "item")
  s$factors[[3]] <- list(name = "length", levels = c("short", "long"),
                         contrasts = list(len = c(-0.5, 0.5)), between = "item")
  s$fixed$coefficients <- list(cond = 0.05, freq = 0.03, len = 0.02)
  d <- simulate_design(s)
  items <- d[!duplicated(d$item), ]
  # The 24 items fall into four cells of six, in item order.
  expect_identical(paste(items$frequency, items$length),
                   rep(c("low short", "low long", "high short", "high long"), each = 6))
})

test_that("joint allocation draws nothing, so it gives the data of one factor of cells", {
  joint <- between_2x2()
  joint$spec_version <- "0.4"
  joint$units$item <- list(n = 3)
  joint$random <- list(subject = list(intercept_sd = 1), item = list(intercept_sd = 0.5))
  # The worked encoding that every version reads: one factor whose levels are the cells, in the
  # same order, carrying both contrast columns.
  cells <- joint
  cells$spec_version <- "0.3"
  cells$factors <- list(list(name = "cell", levels = c("a1.b1", "a1.b2", "a2.b1", "a2.b2"),
                             contrasts = list(a = c(-0.5, -0.5, 0.5, 0.5),
                                              b = c(-0.5, 0.5, -0.5, 0.5)),
                             between = "subject"))
  dj <- simulate_design(joint)
  dc <- simulate_design(cells)
  expect_identical(cells_of(dj, "."), dc$cell)
  expect_identical(dj$y, dc$y)
})

test_that("the parity cases for joint allocation hold the cells they describe", {
  joint <- simulate_design(parity_case("two_between_joint"))
  cells <- simulate_design(parity_case("between_cells_2x2"))
  expect_identical(cells_of(joint, "."), cells$cell)
  expect_identical(joint$score, cells$score)

  d <- simulate_design(parity_case("three_by_two_joint"))
  subjects <- d[!duplicated(d$subject), ]
  expect_identical(paste(subjects$dose, subjects$group),
                   rep(c("d0 g1", "d0 g2", "d1 g1", "d1 g2", "d2 g1", "d2 g2"), each = 2))
  items <- d[d$subject == 1L & d$cond == "x", ]
  expect_identical(paste(items$freq, items$len),
                   rep(c("lo short", "lo long", "hi short", "hi long"), each = 2))
})

test_that("the cells of a joint allocation nest in a grouping factor over the same unit", {
  s <- load_spec(pilotr_example("nested_clusters"))
  s$spec_version <- "0.4"
  s$factors <- between_2x2()$factors
  s$fixed$coefficients <- list(a = 0.3, b = 0.2)
  cells_per_site <- function(s) {
    d <- first_rows(simulate_design(s))
    as.vector(tapply(cells_of(d), d$site, function(x) length(unique(x))))
  }
  # 12 sites, a multiple of the 4 cells: each site lies wholly in one cell.
  expect_identical(cells_per_site(s), rep(1L, 12))
  # 6 sites: sites 2 and 5 straddle two cells each.
  s$random$site$n <- 6
  expect_identical(cells_per_site(s), c(1L, 2L, 1L, 1L, 2L, 1L))
})

test_that("declaring 0.4 changes nothing for one between factor per unit", {
  # With one factor, the joint rule is the block rule of earlier versions, and it draws nothing.
  # tools/parity/golden.json holds every parity case to the same.
  for (ex in c("between_2group_gaussian", "nested_clusters")) {
    s <- load_spec(pilotr_example(ex))
    d <- simulate_design(s)
    s$spec_version <- "0.4"
    expect_identical(simulate_design(s), d)
  }
  s <- parity_case("two_list_counterbalanced")
  d <- simulate_design(s)
  s$spec_version <- "0.4"
  expect_identical(simulate_design(s), d)
})
