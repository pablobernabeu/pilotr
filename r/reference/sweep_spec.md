# Sweep an analysis over one field of a design specification

Vary a single field of a specification across a set of values, run an
analysis at each, and bind the results into one data frame. Sample size
is the axis users sweep most often, and
[`power_curve_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_curve_mixed.md)
and
[`precision_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_curve.md)
are wrappers around this for it, but any field can be swept, including
an effect size, a random-effect standard deviation, a residual standard
deviation, or the number of items per subject.

## Usage

``` r
sweep_spec(
  spec,
  path,
  values,
  fn,
  ...,
  .name = NULL,
  .seeds = c("independent", "shared")
)
```

## Arguments

- spec:

  A design specification (path or list).

- path:

  The field to vary, as `"units$subject$n"` or
  `c("units", "subject", "n")`.

- values:

  A vector or list of values to set the field to, one grid point each.
  [`design_conditions()`](https://pablobernabeu.github.io/pilotr/r/reference/design_conditions.md)
  builds such a list for an effect-size sweep.

- fn:

  The analysis to run at each grid point, for example
  [`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md)
  or
  [`precision_design()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md).
  It is called as `fn(spec, ...)`.

- ...:

  Further arguments passed to `fn` at every grid point. A cluster passed
  as `cl` to an analysis that takes one, such as
  [`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md),
  serves every point, so that one pool of workers runs the whole sweep.

- .name:

  Name for the column recording the swept value. Defaults to the last
  element of `path`.

- .seeds:

  How the grid points are seeded. `"independent"`, the default, gives
  point `i` of `J` the seed `replicate_seeds(spec$seed, J)[i]`, and
  `"shared"` gives every point the specification's own seed, as sweeps
  did up to 0.3.1 (see Details).

## Value

A data frame binding the results, with the swept value as the leading
column, or the grid index when the values are lists.
[`solve_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/solve_curve.md)
reads that leading column by default, so a sweep over single values goes
straight into a solve.

## Details

The specification is validated once, before the sweep, so a mistake in
it is reported before any fitting starts rather than repeated at every
grid point.

A value may be a scalar, which replaces the addressed field, or a list,
which replaces it wholesale. The conditions
[`design_conditions()`](https://pablobernabeu.github.io/pilotr/r/reference/design_conditions.md)
builds are the exception, because each names only the effects it varies.
Each condition is merged into the addressed list with
[`utils::modifyList()`](https://rdrr.io/r/utils/modifyList.html), so the
coefficients it names take its values and every other coefficient keeps
the specification's, in the specification's order. That is the order in
which the linear predictor is summed, so a condition simulates exactly
the data of the specification edited by hand. A condition may name only
coefficients the specification already has, and any other name is
refused.

Each grid point draws its own replicates. Point `i` of a sweep over `J`
values takes `replicate_seeds(spec$seed, J)[i]` as its seed, and the
analysis derives its replicates' seeds from that, so no two points share
their random numbers.
[`solve_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/solve_curve.md)
fits the points of a curve as independent binomial counts, and its
interval keeps its nominal coverage only when they are. The rule is
prefix-stable: appending values to a sweep leaves the earlier points'
results as they were. It depends on order, though, since a point's seed
follows its position in `values`, so reordering the values, or inserting
one, changes the results at every point that moves. A sweep over `seed`
itself sets each point's seed to the value swept.

`.seeds = "shared"` gives every point the specification's own seed, as
every sweep did up to 0.3.1. The points then share their noise, so the
curve looks smoother than its Monte Carlo error allows and an interval
fitted to it is too narrow. Over 1,000 emulated curves for each of
twelve settings on a tight grid of sample sizes, the nominal 95%
interval of
[`solve_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/solve_curve.md)
covered the true size in 87.8% of solves with shared seeds and in 96.5%
with independent ones. On an effect-size axis, where the design and so
every draw stay the same from point to point, the figures were 70.3% and
95.7%. `tools/calibration/solve_curve_coverage.R` in the source
repository produces them.

The result of `fn` is coerced to a data frame: a `pilotr_power` object
becomes one row per focal effect, a data frame is used as it stands, and
a plain list becomes a single row. The swept value is added as a leading
column named by `.name`. A list has no single value to record, so a
sweep over lists, the conditions
[`design_conditions()`](https://pablobernabeu.github.io/pilotr/r/reference/design_conditions.md)
builds included, records the grid index there.
[`solve_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/solve_curve.md)
then needs `x` to name the column that holds the effect. To vary a
single effect, sweep its own coefficient. For example,
`sweep_spec(spec, "fixed$coefficients$cond", c(0, 0.02, 0.04), power_mixed)`
keeps the other coefficients and names the leading column `cond`, so
`solve_curve(curve, target = 0.8, effect = "cond", transform = "identity")`
returns the minimum detectable effect.

## See also

[`design_conditions()`](https://pablobernabeu.github.io/pilotr/r/reference/design_conditions.md)
to build effect-size grids,
[`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md)
and
[`precision_design()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md)
for the analyses usually swept, and
[`solve_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/solve_curve.md)
to solve the resulting curve for the swept value that meets a target.

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

  # Sample size, the same sweep power_curve_mixed() performs.
  sweep_spec(spec, "units$subject$n", c(12, 18), power_mixed, n_sims = 8)

  # Effect size, which the old curve functions could not reach.
  sweep_spec(spec, "fixed$coefficients",
             design_conditions(effect = c(0.03, 0.06)), power_mixed, n_sims = 8)
}
#>   coefficients effect true power power_mcse   power_lo  power_hi n_significant
#> 1            1 effect 0.00 0.000  0.0000000 0.00000000 0.3244076             0
#> 2            2 effect 0.03 0.000  0.0000000 0.00000000 0.3244076             0
#> 3            3 effect 0.06 0.125  0.1169268 0.02241749 0.4708882             1
#>   type_s type_s_lo type_s_hi   type_m type_m_mcse n_attempted n_returned
#> 1     NA        NA        NA       NA          NA           8          8
#> 2     NA        NA        NA       NA          NA           8          8
#> 3      0         0 0.7934507 1.572547          NA           8          8
#>   n_converged n_singular n_warning                          fitter
#> 1           4          4         0 lmerTest::lmer, Satterthwaite t
#> 2           1          7         0 lmerTest::lmer, Satterthwaite t
#> 3           0          8         0 lmerTest::lmer, Satterthwaite t
# }
```
