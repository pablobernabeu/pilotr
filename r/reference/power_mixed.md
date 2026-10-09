# Simulation-based power and design analysis for a mixed-effects design

For each replicate, simulate from the ground-truth specification, fit
the model the specification implies and test each focal fixed effect.
The fitter follows the response family:
[`lmerTest::lmer()`](https://rdrr.io/pkg/lmerTest/man/lmer.html) with
Satterthwaite p-values for the continuous families, and
[`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) with a Wald z
for accuracy and count responses (see Details). Reports power together
with the Type S and Type M errors of Gelman and Carlin (2014). Requires
the `lme4` and `lmerTest` packages.

## Usage

``` r
power_mixed(
  spec,
  focal = NULL,
  formula = NULL,
  prep = NULL,
  n_sims = 100,
  alpha = 0.05,
  workers = 1
)
```

## Arguments

- spec:

  A design specification (path or list).

- focal:

  The fixed effects to test. `NULL`, the default, tests every
  coefficient in the specification and takes the true values from it. A
  character vector names the effects and leaves the true values unknown,
  which suppresses Type S and Type M. A named numeric vector gives both,
  which is how to test against a value other than the one simulated.
  Interaction effects follow the model's column naming, so a
  specification key `a:b` is the focal name `a_b`.

- formula:

  Optional `lme4` formula; if `NULL` it is derived from the
  specification via
  [`model_formula()`](https://pablobernabeu.github.io/pilotr/r/reference/model_formula.md),
  and the response family chooses the fitter. A formula given here is
  fitted as a linear model whatever the family (see Details).

- prep:

  Optional function mapping a simulated data set to the modelling data;
  if `NULL` it is derived via
  [`model_data()`](https://pablobernabeu.github.io/pilotr/r/reference/model_data.md).

- n_sims:

  Number of Monte Carlo replicates. A power estimate carries a Monte
  Carlo standard error of about `sqrt(p * (1 - p) / n_sims)`, and
  `type_s` and `type_m` average over the significant replicates alone,
  so they settle more slowly still. At least 200 replicates are
  advisable for study planning.

- alpha:

  Two-sided significance level.

- workers:

  Number of local worker processes over which to spread the replicates.
  The default of 1 runs serially. Because the replicate seeds are
  derived once from the specification's seed, any worker count returns
  results identical to a serial run. The mixed-model fits dominate the
  cost, so the speed-up is close to linear in the number of cores.

## Value

An object of class `pilotr_power`, a list whose per-run elements are
`n_sims`, `alpha`, `fitter`, `n_attempted`, `n_returned`, `n_converged`,
`n_singular` and `n_warning`, and whose per-effect elements are vectors
named by focal effect: `power`, `power_mcse`, `power_lo`, `power_hi`,
`n_significant`, `true_effect`, `mean_estimate`, `type_s` and `type_m`.
With a single focal effect each of those has length one, so
`result$power` reads as it always has. `fitter` names the fitter and its
test, such as `"lme4::glmer (binomial), Wald z"`. `mean_estimate` and
`type_m` are `NA` for an `ordinal` or `beta` response fitted on its own
scale, and so is the `type_s` of an interaction (see Details).

`power` is the proportion of significant results among the replicates
that returned an estimate for that effect, not among `n_sims`. The
counts report the fit outcomes separately, because a fit can return a
usable estimate while still being boundary-singular or carrying a
convergence warning: `n_returned` counts replicates that yielded a fit,
`n_converged` those that did so with neither a warning nor a singular
fit, `n_singular` those where
[`lme4::isSingular()`](https://rdrr.io/pkg/lme4/man/isSingular.html) was
true, and `n_warning` those with a warning or optimiser convergence
message other than lme4's singular-fit notice, which `n_singular`
already counts. Singular and warning fits are retained in `power`, since
their fixed-effect estimates remain interpretable and discarding them
would bias the result: singularity is not independent of the variance
estimates that produce it. A large `n_singular` means the model being
fitted is richer than the design can support at that sample size, which
is common in crossed designs (Bates et al., 2015; Matuschek et al.,
2017), and is worth reporting alongside the power.

## Details

`power_mixed()` is not a wrapper around an existing power package: it
runs pilotr's own simulation loop over the portable design
specification. It covers territory pioneered by `simr` (Green and
MacLeod, 2016) and `mixedpower` (Kumle, Vo and Draschkow, 2021), to
which it is indebted. pilotr differs in being driven by the portable
cross-language specification, in reporting the Type S and Type M
design-analysis errors alongside power, and in parallelising its
replicates through the `workers` argument.

The analysis model comes from the specification rather than from this
function. Before 0.3 the formula was written into the source as a
maximal crossed structure, so a design declaring uncorrelated slopes, or
no slopes at all, was nonetheless analysed as though it had them, and a
design with more than one factor was refused. The formula now comes from
[`model_formula()`](https://pablobernabeu.github.io/pilotr/r/reference/model_formula.md)
and the data from
[`model_data()`](https://pablobernabeu.github.io/pilotr/r/reference/model_data.md),
so the analysis matches the process that generated the data. Both can
still be given directly, which is what to do when a deliberately
different analysis model is the point, as when checking how a
misspecified model behaves.

The fitter follows the response family, so that the estimates are on the
scale of the specification's coefficients. A `gaussian`, `lognormal`,
`shifted_lognormal` or `exgaussian` response is fitted by
[`lmerTest::lmer()`](https://rdrr.io/pkg/lmerTest/man/lmer.html), the
two lognormal families on the log scale to which
[`model_data()`](https://pablobernabeu.github.io/pilotr/r/reference/model_data.md)
takes them, and each effect is tested with Satterthwaite's degrees of
freedom. A `bernoulli` or `poisson` response is fitted by
[`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) with the
binomial or Poisson family, on the logit or log scale, and each effect
is tested with a Wald z. A model with no random terms, such as that of a
between-subjects design with one row per subject, is fitted by
[`stats::lm()`](https://rdrr.io/r/stats/lm.html) or
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) in the same way. The
result's `fitter` names the fitter and the test. A Wald z treats the
estimate over its standard error as normal and so ignores the
uncertainty in the variance components. Bolker et al. (2009) recommend
it for the fixed effects of a GLMM without overdispersion, which these
models are, but it rejects too often when the clusters are few. In
binary GLMMs with fewer than 30 clusters, a test with as many degrees of
freedom as observations, close to a z test, rejected more often than its
nominal level (Li and Redden, 2015). The power of a `bernoulli` or
`poisson` design with few subjects or items is therefore somewhat
overstated.

pilotr has no frequentist model for an `ordinal` or `beta` response yet.
Such a response is fitted by the linear model on its own scale, while
the specification's coefficients are on the logit scale. The function
therefore warns, and withholds the mean estimate and Type M (`NA`),
since both would compare the two scales. Type S compares signs alone,
and the link keeps the sign of a main effect unless an interaction
reverses that effect between the levels of another factor. It does not
keep the sign of an interaction. A bounded response compresses
differences near the ends of its range, so an interaction that is
positive on the logit scale can be negative on the response scale. Type
S is therefore withheld for an interaction as well. Power is the rate at
which the linear model's test rejects.
[`generate_design_analysis()`](https://pablobernabeu.github.io/pilotr/r/reference/generate_design_analysis.md)
writes a Bayesian design analysis on the link scale, with the model
[`brms_bridge()`](https://pablobernabeu.github.io/pilotr/r/reference/brms_bridge.md)
derives.

A `formula` given here is fitted as written, by
[`lmerTest::lmer()`](https://rdrr.io/pkg/lmerTest/man/lmer.html), or by
[`stats::lm()`](https://rdrr.io/r/stats/lm.html) when it has no random
terms, whatever the family, since only its author knows the scale of its
response. With such a formula, or a `prep` that changes the response,
and `focal = NULL`, the true values still come from the specification
and are on its scale, which need not be the model's. A named numeric
`focal` gives the true values on the model's scale.

Every reported rate carries its Monte Carlo standard error and a Wilson
interval, because a proportion over a finite number of replicates is an
estimate rather than a fact. At the default 100 replicates a power near
0.5 has a standard error of 0.05.

## References

Gelman, A. and Carlin, J. (2014). Beyond power calculations: Assessing
Type S (sign) and Type M (magnitude) errors. *Perspectives on
Psychological Science*, 9(6), 641-651.
[doi:10.1177/1745691614551642](https://doi.org/10.1177/1745691614551642)

Green, P. and MacLeod, C. J. (2016). SIMR: An R package for power
analysis of generalized linear mixed models by simulation. *Methods in
Ecology and Evolution*, 7(4), 493-498.
[doi:10.1111/2041-210x.12504](https://doi.org/10.1111/2041-210x.12504)

Kumle, L., Vo, M. L.-H. and Draschkow, D. (2021). Estimating power in
(generalized) linear mixed models: An open introduction and tutorial in
R. *Behavior Research Methods*, 53, 2528-2543.
[doi:10.3758/s13428-021-01546-0](https://doi.org/10.3758/s13428-021-01546-0)

Bates, D., Kliegl, R., Vasishth, S. and Baayen, H. (2015). Parsimonious
mixed models. *arXiv*.
[doi:10.48550/arXiv.1506.04967](https://doi.org/10.48550/arXiv.1506.04967)

Matuschek, H., Kliegl, R., Vasishth, S., Baayen, H. and Bates, D.
(2017). Balancing Type I error and power in linear mixed models.
*Journal of Memory and Language*, 94, 305-315.
[doi:10.1016/j.jml.2017.01.001](https://doi.org/10.1016/j.jml.2017.01.001)

Bolker, B. M., Brooks, M. E., Clark, C. J., Geange, S. W., Poulsen, J.
R., Stevens, M. H. H. and White, J.-S. S. (2009). Generalized linear
mixed models: A practical guide for ecology and evolution. *Trends in
Ecology & Evolution*, 24(3), 127-135.
[doi:10.1016/j.tree.2008.10.008](https://doi.org/10.1016/j.tree.2008.10.008)

Li, P. and Redden, D. T. (2015). Comparing denominator degrees of
freedom approximations for the generalized linear mixed model in
analyzing binary outcome in small sample cluster-randomized trials. *BMC
Medical Research Methodology*, 15, 38.
[doi:10.1186/s12874-015-0026-x](https://doi.org/10.1186/s12874-015-0026-x)

## See also

[`precision_design()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md)
for the interval-width and ROPE analogue, and
[`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md)
to run this over a grid of sample sizes or effect sizes.

## Examples

``` r
# \donttest{
if (requireNamespace("lme4", quietly = TRUE) &&
    requireNamespace("lmerTest", quietly = TRUE)) {
  spec <- build_spec(list(name = "p", seed = 1, design_kind = "within",
    include_items = TRUE, n_subject = 12, n_item = 12, factor_name = "cond",
    lev1 = "a", lev2 = "b", intercept = 6, effect = 0.05,
    subj_int_sd = 0.12, subj_slope_sd = 0.04, subj_corr = 0.2,
    item_int_sd = 0.08, item_slope_sd = 0.02, item_corr = -0.1,
    family = "shifted_lognormal", resp_name = "", sigma = 0.3, shift = 200))
  # n_sims is small so the example runs quickly. Use 200 or more for real planning.
  power_mixed(spec, n_sims = 10)
}
#> Simulation-based power over 10 replicates (alpha = 0.05)
#>   fitter: lmerTest::lmer, Satterthwaite t
#>   fits: 10 attempted, 10 returned, 0 converged cleanly, 10 singular, 0 with warnings
#>   note: many fits were boundary-singular, so this model is richer than the design supports
#>  effect true power  mcse           ci95 n_sig type_s type_m
#>  effect 0.05   0.2 0.126 [0.057, 0.510]     2      0  1.894
# }
```
