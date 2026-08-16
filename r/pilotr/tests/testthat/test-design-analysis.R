# Tests for the emitted Bayesian design-analysis script. Nothing here fits a model: the whole
# point of generate_design_analysis() is that it templates a script, so what can be checked
# without brms or Stan is that every emitted part is valid R and that the record it writes
# identifies the run that produced it.

# The parts of an array = "slurm" emission, split at the banners. The middle part is bash.
da_parts <- function(script) {
  lines <- strsplit(script, "\n", fixed = TRUE)[[1]]
  at <- grep("^# ===== FILE", lines)
  stopifnot(length(at) == 3L)
  list(analysis = lines[(at[1] + 1L):(at[2] - 1L)],
       slurm = lines[(at[2] + 1L):(at[3] - 1L)],
       aggregate = lines[(at[3] + 1L):length(lines)])
}

# One per-replicate row as the emitted script writes it, so the aggregator can be exercised
# without a sampler.
da_row <- function(rep, rope = 0.05, fingerprint = "0f0f", version = "0.3.0") {
  data.frame(rep = rep, param = "cond", true = 0.05, estimate = 0.052,
             lower = 0.011, upper = 0.093, bf_effect = 12, bf_null = 1 / 12,
             max_rhat = 1.001, p_divergent = 0, verdict = "supported",
             bf_threshold = 10, rope = rope, ci_mass = 0.95,
             max_rhat_limit = 1.01, max_divergent_limit = 0.01,
             spec_fingerprint = fingerprint, pilotr_version = version,
             stringsAsFactors = FALSE)
}

test_that("both emitted R parts are valid R", {
  spec <- load_spec(pilotr_example("crossed_mixed_rt"))
  p <- da_parts(generate_design_analysis(spec, focal = c(cond = 0.05), array = "slurm"))
  expect_type(parse(text = paste(p$analysis, collapse = "\n")), "expression")
  expect_type(parse(text = paste(p$aggregate, collapse = "\n")), "expression")
  expect_true(any(grepl("^#SBATCH --array", p$slurm)))
})

# The wrapper has to be submittable by any user on any SLURM cluster once its two marked
# placeholders are filled in. It previously hard-coded the author's account and cluster
# paths, so no other user's submission could run.
test_that("the slurm wrapper carries placeholders, not the author's cluster", {
  spec <- load_spec(pilotr_example("crossed_mixed_rt"))
  p <- da_parts(generate_design_analysis(spec, focal = c(cond = 0.05), array = "slurm"))
  s <- paste(p$slurm, collapse = "\n")
  expect_false(grepl("educ-intract", s, fixed = TRUE))
  expect_false(grepl("educ1242", s, fixed = TRUE))
  expect_true(grepl("EDIT THESE", s, fixed = TRUE))
  expect_true(any(grepl("^#SBATCH --account=EDIT_ME_ACCOUNT$", p$slurm)))
  expect_true(any(grepl("^PROJECT_DIR=", p$slurm)))
  # The wrapper runs the analysis script saved next to it, not a path fixed at emission.
  expect_true(grepl('"$SCRIPT_DIR/design_analysis.R"', s, fixed = TRUE))
  expect_false(grepl("$HOME/pilotr_toolkit/scripts", s, fixed = TRUE))
})

test_that("a script written to `file` is LF-only with a single trailing newline", {
  spec <- load_spec(pilotr_example("crossed_mixed_rt"))
  f <- tempfile(fileext = ".R"); on.exit(unlink(f))
  out <- generate_design_analysis(spec, focal = c(cond = 0.05), array = "slurm", file = f)
  bytes <- readBin(f, "raw", n = file.size(f) + 16L)
  expect_false(any(bytes == as.raw(0x0D)))            # no CR on any platform
  expect_identical(rawToChar(bytes), out)             # on disk exactly as returned
  expect_true(endsWith(out, "\n"))
  expect_false(endsWith(out, "\n\n"))
})

# The verdict is a pure function of the rule, the gate and the specification, so a record that
# carries only the verdict cannot be told apart from one produced under a different ROPE.
test_that("the emitted record carries the rule, the gate, the spec and the version", {
  spec <- load_spec(pilotr_example("crossed_mixed_rt"))
  script <- generate_design_analysis(spec, focal = c(cond = 0.05), rule = list(rope = 0.02))
  for (col in c("bf_threshold = bf_threshold", "rope = rope", "ci_mass = ci_mass",
                "max_rhat_limit = max_rhat", "max_divergent_limit = max_divergent",
                "spec_fingerprint = spec_fingerprint", "pilotr_version = pilotr_version"))
    expect_true(grepl(col, script, fixed = TRUE))
  # The fingerprint is computed from the embedded specification, not pasted in as a literal.
  expect_true(grepl("tools::md5sum", script, fixed = TRUE))
})

test_that("the aggregator combines one run and refuses to pool two", {
  spec <- load_spec(pilotr_example("crossed_mixed_rt"))
  agg <- paste(da_parts(generate_design_analysis(spec, focal = c(cond = 0.05),
                                                 array = "slurm"))$aggregate, collapse = "\n")
  dir <- file.path(tempdir(), "da-agg")
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  old <- Sys.getenv("PILOTR_OUTDIR", unset = NA)
  on.exit({
    unlink(dir, recursive = TRUE)
    if (is.na(old)) Sys.unsetenv("PILOTR_OUTDIR") else Sys.setenv(PILOTR_OUTDIR = old)
  }, add = TRUE)
  Sys.setenv(PILOTR_OUTDIR = dir)

  saveRDS(da_row(1), file.path(dir, "design_analysis_rep0001.rds"))
  saveRDS(da_row(2), file.path(dir, "design_analysis_rep0002.rds"))
  out <- capture.output(eval(parse(text = agg), envir = new.env()))
  expect_true(any(grepl("combined 2 replicate file", out, fixed = TRUE)))

  # The same two replicates, but the second was run under a different region of practical
  # equivalence. Binding them would report one decision probability over two rules.
  saveRDS(da_row(2, rope = 0.02), file.path(dir, "design_analysis_rep0002.rds"))
  expect_error(eval(parse(text = agg), envir = new.env()), "come from more than one run")

  # An older per-task file, written before the run identity was recorded.
  saveRDS(da_row(2)[, 1:11], file.path(dir, "design_analysis_rep0002.rds"))
  expect_error(eval(parse(text = agg), envir = new.env()), "do not share a column set")
})
