# The lme4 values that pin the Python mixed-model fit.
#
# python/tests/test_power_mixed.py checks that Python's power_mixed() reaches the REML optimum that
# lme4 reaches on the same data. It fits the first replicate of the shipped crossed_mixed_rt
# example and compares the standard error of the condition effect and the REML log-likelihood with
# the values this script prints. Python fits independent by-subject and by-item intercept and slope
# components, so the reference is lme4's uncorrelated model, (1 + cond || subject) +
# (1 + cond || item), fitted to log(RT - shift) as model_data() prepares it.
#
# The data are simulated from the R sources in this checkout, which the parity harness holds
# identical to the Python simulator's. Rerun the script, and update the values pinned in the test,
# whenever a change to the simulator or to the example moves that replicate's data.
#
# Usage: Rscript tools/parity/lme4_reference.R

root <- normalizePath(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "..", ".."), mustWork = FALSE)
if (!nzchar(root) || is.na(root)) root <- normalizePath(".")

if (!requireNamespace("lme4", quietly = TRUE))
  stop("this script needs lme4, which is not installed", call. = FALSE)

# Source the whole package, as tools/parity/run_r.R does, so that the reference comes from the
# simulator in this checkout.
src <- file.path(root, "r", "pilotr", "R")
for (f in sort(list.files(src, pattern = "\\.R$", full.names = TRUE))) source(f)

spec <- load_spec(file.path(root, "spec", "examples", "crossed_mixed_rt.json"))
# Replicate 0 of power_mixed() in either language, which simulates each replicate at its own seed.
replicate <- spec
replicate[["seed"]] <- replicate_seeds(spec[["seed"]], 1)[1]
d <- model_data(spec, simulate_design(replicate))

# REML, as in both power_mixed() functions. calc.derivs = FALSE matches the R twin's fits and
# changes no estimate, since it only skips the convergence checks on the gradient and Hessian.
# lme4's boundary notice is silenced because the last line below reports singularity.
fit <- suppressMessages(
  lme4::lmer(.y ~ cond + (1 + cond || subject) + (1 + cond || item), data = d, REML = TRUE,
             control = lme4::lmerControl(calc.derivs = FALSE)))

cat(sprintf("lme4 %s on crossed_mixed_rt, replicate 0 (seed %s)\n",
            as.character(utils::packageVersion("lme4")),
            format(replicate[["seed"]], scientific = FALSE)))
cat(sprintf("  estimate of cond        %.6f\n", lme4::fixef(fit)[["cond"]]))
cat(sprintf("  standard error of cond  %.6f\n", sqrt(diag(as.matrix(stats::vcov(fit))))[["cond"]]))
cat(sprintf("  REML log-likelihood     %.4f\n", as.numeric(stats::logLik(fit))))
cat(sprintf("  singular                %s\n", lme4::isSingular(fit)))
