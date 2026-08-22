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
sweep_spec(spec, path, values, fn, ..., .name = NULL)
```

## Arguments

- spec:

  A design specification (path or list).

- path:

  The field to vary, as `"units$subject$n"` or
  `c("units", "subject", "n")`.

- values:

  A vector or list of values to set the field to, one grid point each.

- fn:

  The analysis to run at each grid point, for example
  [`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md)
  or
  [`precision_design()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md).
  It is called as `fn(spec, ...)`.

- ...:

  Further arguments passed to `fn` at every grid point.

- .name:

  Name for the column recording the swept value. Defaults to the last
  element of `path`.

## Value

A data frame binding the results, with the swept value as the leading
column. When the swept values are not scalars, that column holds the
grid index instead.
[`solve_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/solve_curve.md)
reads that leading column, so a sweep goes straight into a solve.

## Details

The specification is validated once, before the sweep, so a mistake in
it is reported before any fitting starts rather than repeated at every
grid point.

A value may be a scalar, which replaces the addressed field, or a list,
which replaces it wholesale. Replacing a whole `fixed$coefficients`
object is how an effect-size sweep works, and
[`design_conditions()`](https://pablobernabeu.github.io/pilotr/r/reference/design_conditions.md)
builds those objects, including a common all-zero condition for
examining behaviour under the null.

The result of `fn` is coerced to a data frame: a `pilotr_power` object
becomes one row per focal effect, a data frame is used as it stands, and
a plain list becomes a single row. The swept value is added as a leading
column, named by `.name` where the value is a scalar.

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
             design_conditions(effect = c(0, 0.03, 0.06)), power_mixed, n_sims = 8)
}
#>   coefficients effect true power power_mcse   power_lo  power_hi n_significant
#> 1            1 effect 0.00 0.000  0.0000000 0.00000000 0.3244076             0
#> 2            2 effect 0.00 0.000  0.0000000 0.00000000 0.3244076             0
#> 3            3 effect 0.03 0.125  0.1169268 0.02241749 0.4708882             1
#> 4            4 effect 0.06 0.375  0.1711633 0.13684429 0.6942576             3
#>   type_s   type_m n_attempted n_returned n_converged n_singular n_warning
#> 1     NA       NA           8          8           0          8         8
#> 2     NA       NA           8          8           0          8         8
#> 3      0 2.553562           8          8           1          7         7
#> 4      0 1.655216           8          8           0          8         8
# }
```
