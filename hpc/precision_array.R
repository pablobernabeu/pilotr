#!/usr/bin/env Rscript
# One SLURM array task: a precision and ROPE design analysis at a single sample size.
# Each task is one point of the sweep that precision_curve() runs over `grid`, and it calls
# precision_design() for that point, so its CSV has precision_curve()'s columns and fit counts.
# The task's replicates run on a fork cluster across its cores. Bind the per-task CSVs together
# with rbind() and pass the curve to solve_curve() for the sample size that reaches a target.
#
# Env: TASK_ID (the index into `grid`, from 0), N_SIMS, ROPE (a half-width such as 0.05, or an
# interval such as -0.02,0.05), OUTDIR and SLURM_CPUS_PER_TASK. pilotr and lme4 come from the
# library R_LIBS names, which precision_array.slurm sets to $PROJECT_DIR/Rlib.
#
# The specification, `grid` and `focal` below are the example's. Edit them for another design,
# and keep the wrapper's --array range at one task per value of `grid`.

suppressPackageStartupMessages(library(pilotr))

task   <- as.integer(Sys.getenv("TASK_ID", "0"))
n_sims <- as.integer(Sys.getenv("N_SIMS", "200"))
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
outdir <- Sys.getenv("OUTDIR", ".")
# ROPE is read as comma-separated numbers, and precision_design() refuses anything but a
# non-negative half-width or an increasing pair before any replicate runs. scan() reads a trailing
# comma as a missing bound, which is refused, where strsplit() would drop it.
rope <- scan(text = Sys.getenv("ROPE", "0.05"), sep = ",", quiet = TRUE)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

spec  <- load_spec(pilotr_example("reading_time_continuous"))
grid  <- c(20, 30, 40, 50, 60, 70, 80, 100)        # subjects, one size per array task
focal <- c(SyntaxPC = 0.10, cond_age = 0.02)       # meaningful (outside ROPE) + negligible (inside)
if (is.na(task) || task < 0L || task >= length(grid))
  stop(sprintf("TASK_ID must index the %d sizes in `grid` from 0; got %s",
               length(grid), Sys.getenv("TASK_ID")), call. = FALSE)
N <- grid[task + 1L]
# The seed sweep_spec() gives this point of a sweep over `grid`, so the points of the combined
# curve are independent, as solve_curve() assumes, and a local precision_curve() over the same
# grid reproduces every task's simulated data. precision_design() derives the replicate seeds
# from it with replicate_seeds().
spec$seed <- replicate_seeds(spec$seed, length(grid))[task + 1L]
spec$units$subject$n <- N

# Forked workers inherit the loaded pilotr and lme4, where a socket cluster would start fresh R
# processes and load both again. Forking needs Linux or macOS, and with one core the task runs
# serially on any platform. The arguments go by name: given second, `rope` would bind to
# `formula`.
t0 <- Sys.time()
cl <- if (ncores > 1L) parallel::makeForkCluster(ncores)
result <- tryCatch(
  precision_design(spec, focal = focal, rope = rope, n_sims = n_sims, cl = cl),
  finally = if (!is.null(cl)) parallel::stopCluster(cl))

# The swept size leads, as in precision_curve(), and solve_curve() reads the leading column as
# the swept value. Fitted summaries depend on the lme4 release and the platform, so each row
# records the versions that produced it.
out <- cbind(n_subject = N, result,
             pilotr_version = as.character(utils::packageVersion("pilotr")),
             lme4_version = as.character(utils::packageVersion("lme4")))
file <- file.path(outdir, sprintf("precision_N%03d.csv", N))
utils::write.csv(out, file, row.names = FALSE)
cat(sprintf("wrote %s | N=%d | returned %d, converged %d of %d | %.0fs\n", file, N,
            result$n_returned[1], result$n_converged[1], n_sims,
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
