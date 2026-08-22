# Power curve over sample size for a mixed-effects design

Sweep the number of subjects and compute mixed-effects power at each.
Pass the result to
[`target_n()`](https://pablobernabeu.github.io/pilotr/r/reference/target_n.md)
for the sample size at which power crosses a target, with an interval on
it. A thin wrapper around
[`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md)
over `units$subject$n`, kept because a sample-size curve is the sweep
users want most often. Requires the `lme4` and `lmerTest` packages.

## Usage

``` r
power_curve_mixed(
  spec,
  subject_ns,
  focal = NULL,
  n_sims = 60,
  alpha = 0.05,
  workers = 1
)
```

## Arguments

- spec:

  A design specification (path or list).

- subject_ns:

  A numeric vector of subject counts to evaluate.

- focal:

  The fixed effects to test, as in
  [`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md).

- n_sims:

  Number of Monte Carlo replicates per point. A power estimate carries a
  Monte Carlo standard error of about `sqrt(p * (1 - p) / n_sims)`,
  reported alongside it as `power_mcse`. The default of 60 gives a
  standard error of 0.065 at a power of 0.5, which is too coarse to
  support a claim about a design; raise it to at least 200 for planning.

- alpha:

  Two-sided significance level.

- workers:

  Number of local worker processes over which to spread the replicates
  at each grid point. The default of 1 runs serially, and any worker
  count returns results identical to a serial run.

## Value

A data frame with one row per sample size and focal effect, with columns
`n_subject`, `effect`, `true`, `power`, `power_mcse`, `power_lo`,
`power_hi`, `n_significant`, `type_s`, `type_m`, and the `n_attempted`,
`n_returned`, `n_converged`, `n_singular` and `n_warning` fit counts.
`n_singular` typically falls as the sample size rises, so reading it
down the sweep shows where the model becomes supportable.

## Details

Like
[`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md),
this runs pilotr's own simulation loop over the portable design
specification rather than wrapping an existing package, and differs from
`simr` (Green and MacLeod, 2016) and `mixedpower` (Kumle, Vo and
Draschkow, 2021) in being driven by that specification, in reporting
Type S and Type M errors, and in built-in parallelisation: with
`workers > 1` a single worker pool is created once and reused across all
sample sizes.

For any axis other than sample size, call
[`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md)
directly. Effect size is the axis a design analysis most often needs
after sample size, and
[`design_conditions()`](https://pablobernabeu.github.io/pilotr/r/reference/design_conditions.md)
builds the coefficient overrides for it.

The curve is the input to
[`target_n()`](https://pablobernabeu.github.io/pilotr/r/reference/target_n.md),
which is where the sample size a preregistration quotes should come
from. Reading the crossing off a plot instead judges points whose Monte
Carlo intervals overlap, and reports the answer without the interval
that goes with it.

## References

Green, P. and MacLeod, C. J. (2016). SIMR: An R package for power
analysis of generalized linear mixed models by simulation. *Methods in
Ecology and Evolution*, 7(4), 493-498.
[doi:10.1111/2041-210x.12504](https://doi.org/10.1111/2041-210x.12504)

Kumle, L., Vo, M. L.-H. and Draschkow, D. (2021). Estimating power in
(generalized) linear mixed models: An open introduction and tutorial in
R. *Behavior Research Methods*, 53, 2528-2543.
[doi:10.3758/s13428-021-01546-0](https://doi.org/10.3758/s13428-021-01546-0)

## See also

[`target_n()`](https://pablobernabeu.github.io/pilotr/r/reference/target_n.md)
to solve the returned curve for a sample size,
[`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md)
for any other axis, and
[`design_conditions()`](https://pablobernabeu.github.io/pilotr/r/reference/design_conditions.md)
for effect-size grids.

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
  power_curve_mixed(spec, subject_ns = c(12, 18), n_sims = 8)
}
#>   n_subject effect true power power_mcse   power_lo  power_hi n_significant
#> 1        12 effect 0.05  0.25  0.1530931 0.07147921 0.5907246             2
#> 2        18 effect 0.05  0.50  0.1767767 0.21521606 0.7847839             4
#>   type_s   type_m n_attempted n_returned n_converged n_singular n_warning
#> 1      0 1.893627           8          8           0          8         8
#> 2      0 1.711562           8          8           1          7         7
# }
```
