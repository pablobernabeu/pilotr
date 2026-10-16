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
  workers = 1,
  seeds = c("independent", "shared"),
  alternative = c("two.sided", "greater", "less")
)
```

## Arguments

- spec:

  A design specification (path or list).

- subject_ns:

  A numeric vector of subject counts to evaluate.

- focal:

  The fixed effects to test, or the named sets of coefficients to test
  jointly, as in
  [`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md).

- n_sims:

  Number of Monte Carlo replicates per point. A power estimate carries a
  Monte Carlo standard error of about `sqrt(p * (1 - p) / n_sims)`,
  reported alongside it as `power_mcse`. The default of 60 gives a
  standard error of 0.065 at a power of 0.5, which is too coarse to
  support a claim about a design; raise it to at least 200 for planning.

- alpha:

  Significance level of the test, two-sided unless `alternative` says
  otherwise.

- workers:

  Number of local worker processes over which to spread the replicates
  at each grid point. The default of 1 runs serially, and any worker
  count returns results identical to a serial run. One pool of that many
  processes serves the whole sweep.

- seeds:

  How the sample sizes are seeded: `"independent"`, the default, or
  `"shared"`, as in
  [`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md)'s
  `.seeds`.

- alternative:

  The alternative hypothesis of every test, `"two.sided"`, `"greater"`
  or `"less"`, as in
  [`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md).

## Value

A data frame with one row per sample size and focal effect, with columns
`n_subject`, `effect`, `true`, `power`, `power_mcse`, `power_lo`,
`power_hi`, `n_significant`, `type_s`, `type_s_lo`, `type_s_hi`,
`type_m`, `type_m_mcse`, the `n_attempted`, `n_returned`, `n_converged`,
`n_singular` and `n_warning` fit counts, and `fitter`, as in
[`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md).
`fitter` names the direction of a one-sided test. With sets of
coefficients in `focal`, `effect` holds each set's name, a column
`coefficients` after it lists the set's coefficients, and `fitter` names
the joint test. `n_singular` typically falls as the sample size rises,
so reading it down the sweep shows where the model becomes supportable.

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

Each sample size draws its own replicates, from the seed
[`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md)
gives its grid point, so the points of the curve are independent, as
[`target_n()`](https://pablobernabeu.github.io/pilotr/r/reference/target_n.md)
assumes when it fits them. `seeds = "shared"` runs every point from the
specification's own seed instead, as up to 0.3.1.

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

In a crossed design, the number of items limits what more subjects can
buy. When the effect of a condition varies from item to item, as a
by-item random slope says it does, part of the standard error of the
effect comes from the items sampled. That part shrinks only as items are
added. As subjects are added, the standard error therefore falls towards
a floor set by the number of items. Power at a fixed number of items
levels off below 1 (Westfall, Kenny and Judd, 2014; Judd, Westfall and
Kenny, 2017). The probit curve that
[`target_n()`](https://pablobernabeu.github.io/pilotr/r/reference/target_n.md)
fits rises towards 1 and cannot bend into such a plateau. A curve that
levels off near the target can then be solved far from the true
crossing. One that levels off just below the target never crosses it,
but Monte Carlo noise can lift a simulated power over the target, and
the curve can then be solved as though it crossed. To compare designs
that trade subjects for items, sweep `units$item$n` with
[`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md)
at each of a few subject counts and solve each curve for the number of
items, as the power vignette shows. Power over items levels off in the
same way, at a ceiling set by the by-subject slope and the number of
subjects. The subject counts must be large enough to put that ceiling
well above the target.

## References

Green, P. and MacLeod, C. J. (2016). SIMR: An R package for power
analysis of generalized linear mixed models by simulation. *Methods in
Ecology and Evolution*, 7(4), 493-498.
[doi:10.1111/2041-210x.12504](https://doi.org/10.1111/2041-210x.12504)

Judd, C. M., Westfall, J. and Kenny, D. A. (2017). Experiments with more
than one random factor: Designs, analytic models, and statistical power.
*Annual Review of Psychology*, 68, 601-625.
[doi:10.1146/annurev-psych-122414-033702](https://doi.org/10.1146/annurev-psych-122414-033702)

Kumle, L., Vo, M. L.-H. and Draschkow, D. (2021). Estimating power in
(generalized) linear mixed models: An open introduction and tutorial in
R. *Behavior Research Methods*, 53, 2528-2543.
[doi:10.3758/s13428-021-01546-0](https://doi.org/10.3758/s13428-021-01546-0)

Westfall, J., Kenny, D. A. and Judd, C. M. (2014). Statistical power and
optimal design in experiments in which samples of participants respond
to samples of stimuli. *Journal of Experimental Psychology: General*,
143(5), 2020-2045.
[doi:10.1037/xge0000014](https://doi.org/10.1037/xge0000014)

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
#> 1        12 effect 0.05 0.000  0.0000000 0.00000000 0.3244076             0
#> 2        18 effect 0.05 0.125  0.1169268 0.02241749 0.4708882             1
#>   type_s type_s_lo type_s_hi   type_m type_m_mcse n_attempted n_returned
#> 1     NA        NA        NA       NA          NA           8          8
#> 2      0         0 0.7934507 1.556333          NA           8          8
#>   n_converged n_singular n_warning                          fitter
#> 1           3          5         0 lmerTest::lmer, Satterthwaite t
#> 2           0          8         0 lmerTest::lmer, Satterthwaite t
# }
```
