# Advanced design analysis for a continuous-predictor, interaction-rich, crossed mixed-effects
# design with by-subject random slopes on continuous predictors. The analysis model is
# AUTO-DERIVED from the spec (model_formula / model_data) -- no hand-coded formula. Shows:
# (1) recovery of the full fixed-effect model; (2) a precision-based ROPE design analysis and
# an N-sweep to find the minimum analysable sample size; (3) a brms bridge for the
# confirmatory Bayesian fit.

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
for (f in c("core.R", "simulate.R", "autoformula.R", "precision.R", "brms_bridge.R"))
  source(file.path(here, "..", "R", f))
suppressPackageStartupMessages(library(lme4))
spec <- load_spec(file.path(here, "..", "..", "..", "spec", "examples", "reading_time_continuous.json"))

# ---- (1) Recovery on a large sample, using the auto-derived model ----
big <- spec; big$units$subject$n <- 100; big$units$item$n <- 80
fit <- suppressWarnings(suppressMessages(
  lmer(model_formula(spec), model_data(spec, simulate_design(big)),
       control = lmerControl(calc.derivs = FALSE))))
truth <- c("(Intercept)" = -1.2, SyntaxPC = 0.10, CoherencePC = -0.06, age = -0.08, cond = 0.04,
           SyntaxPC_age = 0.05, CoherencePC_age = 0.03, cond_age = 0.02)
fe <- fixef(fit)
cat("=== Recovery of fixed effects (auto-derived model; N =", nrow(fit@frame), "rows) ===\n")
cat(sprintf("  %-16s %9s %9s\n", "term", "specified", "recovered"))
for (nm in names(truth)) cat(sprintf("  %-16s %9.3f %9.3f\n", nm, truth[nm], fe[nm]))
vc <- as.data.frame(VarCorr(fit))
gs <- function(v) vc$sdcor[vc$grp == "subject" & vc$var1 == v & is.na(vc$var2)]
cat(sprintf("  [random] by-subject SD SyntaxPC %.3f (spec 0.06), CoherencePC %.3f (spec 0.04); residual %.3f (spec 0.25)\n",
            gs("SyntaxPC"), gs("CoherencePC"), vc$sdcor[vc$grp == "Residual"]))

# ---- (2) Precision-based ROPE design analysis + N-sweep (formula auto-derived) ----
focal <- c(SyntaxPC = 0.10, cond_age = 0.02)   # meaningful (outside ROPE) and negligible (inside)
cat("\n=== Precision-based design analysis: sweep number of subjects (ROPE |beta| < 0.05) ===\n")
curve <- precision_curve(spec, focal, subject_ns = c(20, 40, 60), rope = 0.05, n_sims = 30)
print(curve, digits = 3)
hit <- curve[curve$param == "SyntaxPC" & curve$p_meaningful >= 0.80, "n_subject"]
cat(if (length(hit)) sprintf("  -> SyntaxPC (true 0.10) reaches P(meaningful) >= 0.80 by %d subjects.\n", min(hit))
    else "  -> extend the grid: P(meaningful) >= 0.80 for SyntaxPC not yet reached.\n")

# ---- (3) Bridge to the confirmatory Bayesian fit ----
cat("\n=== brms bridge (confirmatory Bayesian model from the same spec) ===\n")
invisible(brms_bridge(spec))
