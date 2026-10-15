# The SLURM array job in the repository's hpc/ folder. It ran only under the author's account and
# directory layout: it charged a fixed account, wrote its logs and results under the author's
# project paths, ran a copy of its runner kept outside the repository and needed a bootstrap
# script the repository never held. The runner sourced the package files and repeated the
# replicate loop of precision_design(), counting every fit that returned as converged, which is
# the definition 0.3.0 replaced, in columns that precision_curve() does not share. hpc/ sits
# beside the package, not in it, so these tests skip wherever the repository is absent, R CMD
# check included. python/tests/test_hpc.py holds the same guard in the Python suite, which CI runs
# on the whole checkout.

hpc_file <- function(name) {
  path <- testthat::test_path("..", "..", "..", "..", "hpc", name)
  testthat::skip_if_not(file.exists(path), "hpc/ is in the repository, not the package")
  path
}

hpc_lines <- function(name) readLines(hpc_file(name), encoding = "UTF-8", warn = FALSE)

# Every call to `name` anywhere in an expression. An empty argument, such as the row index of
# x[, 1], is never bound to a variable, since evaluating it would fail.
calls_to <- function(e, name) {
  out <- list()
  if (is.call(e)) {
    if (identical(e[[1]], as.name(name))) out <- list(e)
    for (i in seq_along(e)[-1]) if (is.call(e[[i]])) out <- c(out, calls_to(e[[i]], name))
  }
  out
}

test_that("the HPC job names no account, path or layout of the author's cluster", {
  for (name in c("precision_array.slurm", "precision_array.R")) {
    s <- paste(hpc_lines(name), collapse = "\n")
    for (banned in c("educ-intract", "educ1242", "/data/educ", "pilotr_toolkit", "bootstrap_rlib"))
      expect_false(grepl(banned, s, fixed = TRUE), info = paste(name, "holds", banned))
  }
})

# The same placeholders as the wrapper generate_design_analysis() emits, which
# test-design-analysis.R checks.
test_that("the SLURM wrapper carries marked placeholders and runs the runner beside it", {
  lines <- hpc_lines("precision_array.slurm")
  s <- paste(lines, collapse = "\n")
  expect_true(grepl("EDIT THESE", s, fixed = TRUE))
  expect_true(any(grepl("^#SBATCH --account=EDIT_ME_ACCOUNT$", lines)))
  expect_true(any(grepl("^PROJECT_DIR=.*EDIT THIS", lines)))
  # A log path is resolved where the job runs, so an absolute one names somebody's directory.
  expect_false(any(grepl("^#SBATCH +--(output|error)=/", lines)))
  expect_false(any(grepl("^module load R/", lines)))
  expect_true(any(grepl("^SCRIPT_DIR=\\$\\{SLURM_SUBMIT_DIR", lines)))
  expect_true(grepl('Rscript "$SCRIPT_DIR/precision_array.R"', s, fixed = TRUE))
})

test_that("the runner calls precision_design() from the installed package, by name", {
  exprs <- parse(hpc_file("precision_array.R"), keep.source = FALSE)
  attached <- vapply(unlist(lapply(exprs, calls_to, "library"), recursive = FALSE),
                     function(e) as.character(e[[2]]), character(1))
  expect_true("pilotr" %in% attached)
  used <- all.names(exprs)
  expect_false(any(c("source", "sys.source", "mclapply") %in% used))
  found <- unlist(lapply(exprs, calls_to, "precision_design"), recursive = FALSE)
  expect_length(found, 1L)
  # Only the specification goes by position. Given second, a region of practical equivalence
  # would bind to `formula`.
  args <- names(as.list(found[[1]]))[-(1:2)]
  expect_true(all(nzchar(args)))
  expect_true(all(c("focal", "rope", "n_sims", "cl") %in% args))
})

# One array task is one point of the sweep precision_curve() runs over the job's grid, seeded as
# sweep_spec() seeds that point. Its CSV holds that point's row, the swept size leading, with the
# versions of pilotr and lme4 that produced it. The first two seeds of the eight-point grid are
# the seeds of a two-point sweep, so a two-point curve gives the second task's row.
test_that("a task writes the row a local precision_curve() gives its sample size", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  runner <- hpc_file("precision_array.R")
  out <- tempfile("hpc-")
  dir.create(out)
  vars <- c(TASK_ID = "1", N_SIMS = "2", ROPE = "0.05", OUTDIR = out, SLURM_CPUS_PER_TASK = "1")
  old <- Sys.getenv(names(vars), unset = NA, names = TRUE)
  on.exit({
    unlink(out, recursive = TRUE)
    for (v in names(old))
      if (is.na(old[[v]])) Sys.unsetenv(v) else do.call(Sys.setenv, as.list(old[v]))
  }, add = TRUE)
  do.call(Sys.setenv, as.list(vars))

  # Evaluated over this pilotr's namespace, so that the runner's calls reach the code under test
  # whichever pilotr its library() call finds.
  env <- new.env(parent = asNamespace("pilotr"))
  utils::capture.output(sys.source(runner, envir = env))
  written <- utils::read.csv(file.path(out, "precision_N030.csv"), stringsAsFactors = FALSE)

  spec <- load_spec(pilotr_example("reading_time_continuous"))
  local <- precision_curve(spec, focal = c(SyntaxPC = 0.10, cond_age = 0.02),
                           subject_ns = c(20, 30), rope = 0.05, n_sims = 2)
  local <- local[local$n_subject == 30, , drop = FALSE]
  row.names(local) <- NULL
  expect_identical(names(written), c(names(local), "pilotr_version", "lme4_version"))
  # write.csv() keeps 15 significant digits, so the numbers agree to that precision.
  expect_equal(written[names(local)], local)
  expect_identical(unique(written$pilotr_version), as.character(utils::packageVersion("pilotr")))
  expect_identical(unique(written$lme4_version), as.character(utils::packageVersion("lme4")))
})
