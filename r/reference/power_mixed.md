# Simulation-based power and design analysis for a mixed-effects design

For each replicate, simulate from the ground-truth specification, fit
the model the specification implies with `lmerTest`, and test each focal
fixed effect using Satterthwaite p-values. Reports power together with
the Type S and Type M errors of Gelman and Carlin (2014). Requires the
`lme4` and `lmerTest` packages.

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
  [`model_formula()`](https://pablobernabeu.github.io/pilotr/r/reference/model_formula.md).

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
`n_sims`, `alpha`, `n_attempted`, `n_returned`, `n_converged`,
`n_singular` and `n_warning`, and whose per-effect elements are vectors
named by focal effect: `power`, `power_mcse`, `power_lo`, `power_hi`,
`n_significant`, `true_effect`, `mean_estimate`, `type_s` and `type_m`.
With a single focal effect each of those has length one, so
`result$power` reads as it always has.

`power` is the proportion of significant results among the replicates
that returned an estimate for that effect, not among `n_sims`. The
counts report the fit outcomes separately, because a fit can return a
usable estimate while still being boundary-singular or carrying a
convergence warning: `n_returned` counts replicates that yielded a fit,
`n_converged` those that did so with neither a warning nor a singular
fit, `n_singular` those where
[`lme4::isSingular()`](https://rdrr.io/pkg/lme4/man/isSingular.html) was
true, and `n_warning` those with a warning or optimiser convergence
message. Singular and warning fits are retained in `power`, since their
fixed-effect estimates remain interpretable and discarding them would
bias the result: singularity is not independent of the variance
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
#>   fits: 10 attempted, 10 returned, 0 converged cleanly, 10 singular, 10 with warnings
#>   note: many fits were boundary-singular, so this model is richer than the design supports
#>  effect true power  mcse           ci95 n_sig type_s type_m
#>  effect 0.05   0.2 0.126 [0.057, 0.510]     2      0  1.894
# }
```
