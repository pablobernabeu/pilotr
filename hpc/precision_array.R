#!/usr/bin/env Rscript
# One ARC array task: a high-precision ROPE design analysis at a single sample size.
# Embarrassingly parallel -- the task's n_sims replicates run across cpus-per-task cores via
# mclapply. Each replicate simulates from the spec, fits the auto-derived lmer model, and
# records the focal 95% CIs; the task aggregates ROPE decision probabilities and writes one
# CSV per sample size to the project results dir. Combine the CSVs into a full precision curve.
#
# Env: TOOLKIT (toolkit dir), TASK_ID (grid index), N_SIMS, ROPE, OUTDIR, SLURM_CPUS_PER_TASK.

suppressWarnings(suppressMessages({
  TK <- path.expand(Sys.getenv("TOOLKIT", "~/simdgp_toolkit/toolkit"))
  for (f in c("core.R", "simulate.R", "autoformula.R")) source(file.path(TK, "r/simdgp/R", f))
  library(lme4); library(parallel)
}))

task   <- as.integer(Sys.getenv("TASK_ID", "0"))
n_sims <- as.integer(Sys.getenv("N_SIMS", "200"))
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
rope   <- as.numeric(Sys.getenv("ROPE", "0.05"))
outdir <- Sys.getenv("OUTDIR", ".")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

spec <- load_spec(file.path(TK, "spec/examples/reading_time_continuous.json"))
grid <- c(20, 30, 40, 50, 60, 70, 80, 100)        # subjects; one per array task
N <- grid[task + 1L]
spec$units$subject$n <- N
form  <- model_formula(spec)
focal <- c(SyntaxPC = 0.10, cond_age = 0.02)       # meaningful (outside ROPE) + negligible (inside)
base  <- spec$seed

one_sim <- function(i) {
  s <- spec; s$seed <- base + i
  d <- model_data(spec, simulate_design(s))
  fit <- tryCatch(lme4::lmer(form, d, control = lme4::lmerControl(calc.derivs = FALSE)),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  est <- lme4::fixef(fit); se <- sqrt(diag(as.matrix(stats::vcov(fit))))
  vapply(names(focal), function(f) c(est = est[[f]], se = se[[f]]), numeric(2))
}

t0 <- Sys.time()
res <- Filter(Negate(is.null), mclapply(seq_len(n_sims), one_sim, mc.cores = ncores))
nc <- length(res)
agg <- do.call(rbind, lapply(names(focal), function(f) {
  est <- vapply(res, function(r) r["est", f], numeric(1))
  se  <- vapply(res, function(r) r["se",  f], numeric(1))
  lo <- est - 1.96 * se; hi <- est + 1.96 * se
  data.frame(n_subject = N, param = f, true = focal[[f]],
             mean_ci_width = mean(hi - lo),
             p_meaningful = mean(lo > rope | hi < -rope),
             p_equivalent = mean(lo > -rope & hi < rope),
             n_converged = nc, n_sims = n_sims)
}))
out <- file.path(outdir, sprintf("precision_N%03d.csv", N))
write.csv(agg, out, row.names = FALSE)
cat(sprintf("wrote %s | N=%d | converged %d/%d | %.0fs\n", out, N, nc, n_sims,
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
