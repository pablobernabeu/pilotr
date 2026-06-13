# Advanced design analysis for a continuous-predictor, interaction-rich, crossed mixed-effects
# design with by-subject random slopes on continuous predictors. Demonstrates: (1) recovery of
# the full fixed-effect model; (2) a precision-based ROPE design analysis and an N-sweep to
# find the minimum analysable sample size; (3) a brms bridge for the confirmatory Bayesian fit.

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
for (f in c("core.R", "simulate.R", "precision.R", "brms_bridge.R")) source(file.path(here, "..", "R", f))
suppressPackageStartupMessages(library(lme4))
SPEC <- file.path(here, "..", "..", "..", "spec", "examples", "reading_time_continuous.json")

# interactions entered as explicit product columns so coefficient names are exact
FORM <- logRT ~ SyntaxPC + CoherencePC + age + cond + Syntax_age + Coher_age + narr_age +
  (1 + SyntaxPC + CoherencePC | subject) + (1 | item)
prep <- function(d) {
  d$cond <- ifelse(d$narration == "on", 0.5, -0.5)
  d$logRT <- log(d$reading_time_per_word)
  d$Syntax_age <- d$SyntaxPC * d$age; d$Coher_age <- d$CoherencePC * d$age; d$narr_age <- d$cond * d$age
  d
}

# ---- (1) Recovery on a large sample ----
spec <- load_spec(SPEC)
big <- spec; big$units$subject$n <- 100; big$units$item$n <- 80
fit <- suppressWarnings(suppressMessages(lmer(FORM, prep(simulate_design(big)),
                                              control = lmerControl(calc.derivs = FALSE))))
truth <- c("(Intercept)" = -1.2, SyntaxPC = 0.10, CoherencePC = -0.06, age = -0.08, cond = 0.04,
           Syntax_age = 0.05, Coher_age = 0.03, narr_age = 0.02)
fe <- fixef(fit)
cat("=== Recovery of fixed effects (N =", nrow(fit@frame), "rows) ===\n")
cat(sprintf("  %-14s %9s %9s\n", "term", "specified", "recovered"))
for (nm in names(truth)) cat(sprintf("  %-14s %9.3f %9.3f\n", nm, truth[nm], fe[nm]))
vc <- as.data.frame(VarCorr(fit))
gs <- function(v) vc$sdcor[vc$grp == "subject" & vc$var1 == v & is.na(vc$var2)]
cat(sprintf("  [random] by-subject SD SyntaxPC %.3f (spec 0.06), CoherencePC %.3f (spec 0.04); residual %.3f (spec 0.25)\n",
            gs("SyntaxPC"), gs("CoherencePC"), vc$sdcor[vc$grp == "Residual"]))

# ---- (2) Precision-based ROPE design analysis + N-sweep ----
focal <- c(SyntaxPC = 0.10, narr_age = 0.02)   # meaningful (outside ROPE) and negligible (inside)
cat("\n=== Precision-based design analysis: sweep number of subjects (ROPE |beta| < 0.05) ===\n")
curve <- precision_curve(spec, FORM, focal, subject_ns = c(20, 40, 60), rope = 0.05, n_sims = 30, prep = prep)
print(curve, digits = 3)
meaningful <- curve[curve$param == "SyntaxPC" & curve$p_meaningful >= 0.80, "n_subject"]
cat(if (length(meaningful))
      sprintf("  -> SyntaxPC (true 0.10) reaches P(meaningful) >= 0.80 by %d subjects.\n", min(meaningful))
    else "  -> extend the grid: P(meaningful) >= 0.80 for SyntaxPC not yet reached.\n")

# ---- (3) Bridge to the confirmatory Bayesian fit ----
cat("\n=== brms bridge (confirmatory Bayesian model from the same spec) ===\n")
invisible(brms_bridge(spec))
