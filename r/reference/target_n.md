# Solve a power curve for the sample size that reaches a target power

The sample-size case of
[`solve_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/solve_curve.md),
and the number a power analysis is usually run to obtain. Takes the
curve a sweep over sample size has produced and returns the size at
which power reaches `target`, rounded up to a whole number of units
alongside the exact solution.

## Usage

``` r
target_n(curve, target = 0.8, ...)
```

## Arguments

- curve:

  A power curve, as returned by
  [`power_curve_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_curve_mixed.md)
  or by
  [`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md)
  over `units$subject$n`.

- target:

  The power to reach. Defaults to 0.8, the convention this package's
  plots draw a line at.

- ...:

  Further arguments passed to
  [`solve_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/solve_curve.md),
  such as `effect` to pick one focal effect out of a curve holding
  several, or `level` for the interval.

## Value

The list
[`solve_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/solve_curve.md)
returns, with `n`, `n_lo` and `n_hi` added: `value`, `lo` and `hi`
rounded up to whole numbers.

## Details

Everything
[`solve_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/solve_curve.md)
does applies here, including its refusals: a curve that never reaches
the target within the sizes it swept is refused outright, and the
reported interval can extend past the largest size simulated, which
means the sweep was too narrow to settle the question.

The whole-number fields round up rather than to nearest, because a
design cannot recruit a fraction of a subject and rounding down would
leave the study short of the target it was sized for.

## See also

[`solve_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/solve_curve.md),
which this wraps, and
[`power_curve_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_curve_mixed.md)
for the curve.

## Examples

``` r
spec <- build_spec(list(name = "s", seed = 1, design_kind = "between", n_subject = 40,
  factor_name = "group", lev1 = "a", lev2 = "b", intercept = 0, effect = 0.7,
  family = "gaussian", resp_name = "score", sigma = 1))
# n_sims is small so the example runs quickly. Use 200 or more for real planning.
curve <- sweep_spec(spec, "units$subject$n", c(20, 40, 60, 80), power_design, n_sims = 50)
solved <- target_n(curve)
unlist(solved[c("n", "n_lo", "n_hi")])
#>    n n_lo n_hi 
#>   63   53   74 
```
