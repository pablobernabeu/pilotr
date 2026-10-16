# Precision and ROPE curve over sample size

Sweep the number of subjects and report the ROPE decision probabilities
at each size. Pass the result to
[`solve_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/solve_curve.md)
for the minimum analysable *N* at which a focal effect reaches a
determinate decision with a target probability, such as 0.90, together
with an interval on it. Calls
[`precision_design()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md)
and so requires the `lme4` package.

## Usage

``` r
precision_curve(
  spec,
  focal = NULL,
  subject_ns,
  formula = NULL,
  prep = NULL,
  rope = 0.05,
  n_sims = 60,
  workers = 1,
  seeds = c("independent", "shared"),
  fallback = c("singular", "warning")
)
```

## Arguments

- spec:

  A design specification (path or list).

- focal:

  The focal effects, as in
  [`precision_design()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md).
  `NULL` uses every coefficient in the specification.

- subject_ns:

  A numeric vector of subject counts to evaluate.

- formula:

  Optional `lme4` formula; if `NULL` it is derived via
  [`model_formula()`](https://pablobernabeu.github.io/pilotr/r/reference/model_formula.md),
  and the response family chooses the fitter, as in
  [`precision_design()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md).
  A ladder of models, as in
  [`precision_design()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md),
  simulates a simplification rule at every sample size.

- prep:

  Optional data-preparation function; if `NULL` it is derived via
  [`model_data()`](https://pablobernabeu.github.io/pilotr/r/reference/model_data.md).

- rope:

  The region of practical equivalence, as in
  [`precision_design()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md):
  a non-negative half-width `r` for the region from `-r` to `r`, or a
  pair `c(lower, upper)` with `lower < upper`. Keep the effects worth
  detecting clearly outside the region, because the probability of a
  determinate meaningful decision about an effect inside it cannot rise
  above 0.5 however large the sample, so the curve would fall with *N*
  rather than rise.

- n_sims:

  Number of Monte Carlo replicates per sample size. `p_meaningful` and
  `p_equivalent` are proportions over the replicates that produced an
  estimate, so each is reported with its Monte Carlo standard error and
  Wilson interval. The default of 60 gives a standard error of 0.065 at
  a rate of 0.5, which is too coarse to support a claim about a design;
  raise it to at least 200 for planning.

- workers:

  Number of local worker processes over which to spread the replicates
  at each sample size. The default of 1 runs serially, and any worker
  count returns results identical to a serial run. One pool of that many
  processes serves the whole sweep.

- seeds:

  How the sample sizes are seeded: `"independent"`, the default, or
  `"shared"`, as in
  [`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md)'s
  `.seeds`.

- fallback:

  The rule of a ladder, `"singular"` or `"warning"`, as in
  [`precision_design()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md).

## Value

A data frame with one row per focal effect and sample size, adding an
`n_subject` column to the columns returned by
[`precision_design()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md),
including the Monte Carlo standard errors, the Wilson interval bounds,
and the `n_returned`, `n_converged`, `n_singular` and `n_warning` fit
counts, and for a ladder the share of each rung at each sample size.

## Details

A thin wrapper around
[`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md)
over `units$subject$n`, kept because a sample-size curve is the sweep
users want most often. For any other axis, call
[`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md)
directly; effect size is the axis a design analysis most often needs
next, and
[`design_conditions()`](https://pablobernabeu.github.io/pilotr/r/reference/design_conditions.md)
builds the coefficient overrides for it.

Each sample size draws its own replicates, from the seed
[`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md)
gives its grid point, so the points of the curve are independent, as
[`solve_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/solve_curve.md)
assumes when it fits them. `seeds = "shared"` runs every point from the
specification's own seed instead, as up to 0.3.1. With `workers > 1`,
one worker pool is created and reused across all sample sizes.

Reading the crossing off the returned curve, or off a plot of it, is
what
[`solve_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/solve_curve.md)
replaces. At the replicate counts these runs are usually given,
neighbouring points on the curve are not significantly different from
one another, so a sample size judged by eye is a point estimate with an
unstated and often wide uncertainty behind it.

## See also

[`solve_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/solve_curve.md)
to solve the returned curve for a target decision probability,
[`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md)
for any other axis, and
[`design_conditions()`](https://pablobernabeu.github.io/pilotr/r/reference/design_conditions.md)
for effect-size grids.

## Examples

``` r
# \donttest{
if (requireNamespace("lme4", quietly = TRUE)) {
  spec <- build_spec(list(name = "pr", seed = 1, design_kind = "within",
    include_items = TRUE, n_subject = 12, n_item = 12, factor_name = "cond",
    lev1 = "a", lev2 = "b", intercept = 6, effect = 0.05,
    subj_int_sd = 0.12, subj_slope_sd = 0.04, subj_corr = 0.2,
    item_int_sd = 0.08, item_slope_sd = 0.02, item_corr = -0.1,
    family = "shifted_lognormal", resp_name = "", sigma = 0.3, shift = 200))
  # n_sims is small so the example runs quickly. Use 200 or more for real planning.
  precision_curve(spec, focal = c(effect = 0.05), subject_ns = c(12, 18), rope = 0.02,
                  n_sims = 8)
}
#>   n_subject  param true mean_ci_width p_meaningful p_meaningful_mcse
#> 1        12 effect 0.05     0.1901667            0                 0
#> 2        18 effect 0.05     0.1308207            0                 0
#>   p_meaningful_lo p_meaningful_hi p_equivalent p_equivalent_mcse
#> 1               0       0.3244076            0                 0
#> 2               0       0.3244076            0                 0
#>   p_equivalent_lo p_equivalent_hi n_attempted n_returned n_converged n_singular
#> 1               0       0.3244076           8          8           3          5
#> 2               0       0.3244076           8          8           0          8
#>   n_warning             fitter
#> 1         0 lme4::lmer, Wald z
#> 2         0 lme4::lmer, Wald z
# }
```
