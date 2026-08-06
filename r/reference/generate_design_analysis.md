# Generate a Bayesian design-analysis script from a specification

Emit a runnable R script that simulates from a design specification,
fits the confirmatory Bayesian model that
[`brms_bridge()`](https://pablobernabeu.github.io/pilotr/r/reference/brms_bridge.md)
derives for it, and decides about each focal effect on a Savage-Dickey
Bayes factor together with a highest-density interval tested against a
region of practical equivalence. The script is returned as a character
string and nothing is fitted here, so the function needs neither `brms`
nor Stan installed.

## Usage

``` r
generate_design_analysis(
  spec,
  focal,
  rule = list(bf = 10, rope = 0.05, ci = 0.95),
  engine = "brms",
  gate = list(max_rhat = 1.01, max_divergent = 0.01),
  array = c("none", "slurm"),
  file = NULL
)
```

## Arguments

- spec:

  A design specification (path or list).

- focal:

  A character vector of focal coefficient names, or a named numeric
  vector mapping those names to their true values. The names are the
  coefficients of the emitted model, so an interaction is written `a:b`
  as in the specification, not `a_b` as in the `lme4` formula from
  [`model_formula()`](https://pablobernabeu.github.io/pilotr/r/reference/model_formula.md).
  A name that is not a fixed coefficient of the design is reported as a
  warning here rather than left to fail inside the script.

- rule:

  The decision rule, a named list with elements `bf` (the Bayes-factor
  threshold, applied in both directions), `rope` (the half-width of the
  region of practical equivalence, on the scale of the model's
  coefficients) and `ci` (the mass of the highest-density interval). A
  missing element takes its default. Setting `rope` to 0 leaves a
  verdict of `"null"` unreachable, since no interval of positive width
  lies inside a region of no width, which turns the rule into a Bayes
  factor plus a sign requirement.

- engine:

  The fitting engine for the emitted script. Only `"brms"` is supported.

- gate:

  The convergence gate, a named list with elements `max_rhat` (the
  largest acceptable R-hat over all parameters) and `max_divergent` (the
  largest acceptable share of post-warmup draws ending in a divergent
  transition). A missing element takes its default.

- array:

  Either `"none"` for the analysis script alone, or `"slurm"` to append
  a SLURM array wrapper that runs one replicate per task and an
  aggregator that combines the per-task results. The three parts are
  separated by `# ===== FILE n of 3` banners and are meant to be split
  into three files, since the middle part is shell and so is not valid
  R.

- file:

  Optional path. When given, the script is also written there with
  [`writeLines()`](https://rdrr.io/r/base/writeLines.html) and returned
  invisibly.

## Value

A length-one character string holding the emitted script, invisibly when
`file` is given. For `array = "slurm"` the string holds three
banner-separated parts, an R analysis script, a bash array wrapper, and
an R aggregator.

## Details

The function emits a script rather than running the analysis, for two
reasons that both come down to where a Stan model can be built. The
no-code application ships as a webR build that runs entirely in the
browser, and Stan compiles C++ at fit time, so a Bayesian fit cannot run
there at all. A function that produced its analysis only when `brms` was
present would then be missing from the build most users meet first. The
second reason is maintenance. Depending on `brms` would pull a Stan
toolchain into this package's own test and check matrix, and one
maintainer cannot keep that working across the platforms CRAN builds on.
Emitting a script keeps the analysis reproducible and open to inspection
while leaving the fit where a compiler is available.

The verdict rests on two criteria that answer different questions. A
Bayes factor compares the null with the alternative and so reports which
of the two the data favour (Kass and Raftery, 1995), computed here as
the Savage-Dickey density ratio, the ratio of prior to posterior density
at zero (Wagenmakers et al., 2010). The interval against a region of
practical equivalence asks instead whether the effect is large enough to
matter (Kruschke, 2018). Requiring the two to agree makes `"supported"`
and `"null"` harder to reach than either criterion alone would, and it
leaves the third answer of `"inconclusive"` available when they
disagree, which is the answer a small pilot most often deserves.

Because the Savage-Dickey ratio is read off the prior as well as the
posterior, the emitted `brm()` call sets `sample_prior = "yes"`, and the
Bayes factor it produces is a statement about the prior that
[`brms_bridge()`](https://pablobernabeu.github.io/pilotr/r/reference/brms_bridge.md)
supplies as much as about the data. Widening that prior moves the factor
towards the null.

The convergence gate is checked before any verdict is formed. R-hat is
the rank-normalised version of Vehtari et al. (2021), for which a limit
near 1.01 is appropriate rather than the older 1.1, and divergent
transitions are counted as a share of post-warmup draws so that the
threshold means the same thing whatever the run length. When either
limit is exceeded the emitted script reports `NA` for every verdict and
prints why, since a conclusion from a fit that has not converged carries
the authority of a number without the sampling behind it.

## References

Kass, R. E. and Raftery, A. E. (1995). Bayes factors. *Journal of the
American Statistical Association*, 90(430), 773-795.
[doi:10.1080/01621459.1995.10476572](https://doi.org/10.1080/01621459.1995.10476572)

Kruschke, J. K. (2018). Rejecting or accepting parameter values in
Bayesian estimation. *Advances in Methods and Practices in Psychological
Science*, 1(2), 270-280.
[doi:10.1177/2515245918771304](https://doi.org/10.1177/2515245918771304)

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and Burkner, P.-C.
(2021). Rank-normalization, folding, and localization: An improved R-hat
for assessing convergence of MCMC. *Bayesian Analysis*, 16(2), 667-718.
[doi:10.1214/20-BA1221](https://doi.org/10.1214/20-BA1221)

Wagenmakers, E.-J., Lodewyckx, T., Kuriyal, H. and Grasman, R. (2010).
Bayesian hypothesis testing for psychologists: A tutorial on the
Savage-Dickey method. *Cognitive Psychology*, 60(3), 158-189.
[doi:10.1016/j.cogpsych.2009.12.001](https://doi.org/10.1016/j.cogpsych.2009.12.001)

## See also

[`brms_bridge()`](https://pablobernabeu.github.io/pilotr/r/reference/brms_bridge.md)
for the model the script fits,
[`precision_design()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md)
for the frequentist analogue of the interval criterion, which runs in
place, and
[`generate_r_script()`](https://pablobernabeu.github.io/pilotr/r/reference/generate_r_script.md)
for the simulation-only script.

## Examples

``` r
# Emitting the script needs neither brms nor Stan, which is what lets it work in the
# browser build of the no-code app.
spec   <- load_spec(pilotr_example("crossed_mixed_rt"))
script <- generate_design_analysis(spec, focal = c(cond = 0.05))
cat(head(strsplit(script, "\n")[[1]], 15), sep = "\n")
#> #!/usr/bin/env Rscript
#> # ---------------------------------------------------------------------------
#> # Bayesian design analysis for the pilotr design 'crossed_mixed_rt'.
#> #
#> # One run is one replicate. The script simulates from the specification embedded
#> # below, fits the confirmatory model with brms, and decides about each focal
#> # effect on a Savage-Dickey Bayes factor together with a highest-density interval
#> # against a region of practical equivalence. No verdict is reported unless the
#> # sampler has converged.
#> #
#> # Set PILOTR_REP to a replicate index to reseed the design, and PILOTR_OUTDIR to a
#> # directory to have the replicate's summary written there as an RDS. Leave both
#> # unset for a single interactive run.
#> # ---------------------------------------------------------------------------
#> 

# A stricter rule, with a wrapper for a SLURM array and its aggregator.
cluster <- generate_design_analysis(spec, focal = c(cond = 0.05),
                                    rule = list(bf = 30, rope = 0.02),
                                    array = "slurm")
cat(grep("^# ===== FILE", strsplit(cluster, "\n")[[1]], value = TRUE), sep = "\n")
#> # ===== FILE 1 of 3: design_analysis.R (R) ====================================
#> # ===== FILE 2 of 3: design_analysis.slurm (bash) =============================
#> # ===== FILE 3 of 3: aggregate_design_analysis.R (R) ==========================
```
