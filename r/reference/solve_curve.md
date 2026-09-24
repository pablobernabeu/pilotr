# Solve a simulated design curve for the value that meets a target

Take the curve a sweep has already produced, fit the decision rate
against the swept value, and solve for the value at which the rate meets
a target. The solved value comes with a confidence interval, because a
point read off a simulated curve without one repeats the overconfidence
that design analysis exists to expose.

## Usage

``` r
solve_curve(
  curve,
  target,
  x = NULL,
  y = NULL,
  n = NULL,
  effect = NULL,
  transform = "sqrt",
  level = 0.95
)
```

## Arguments

- curve:

  A curve, as returned by
  [`power_curve_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_curve_mixed.md),
  [`precision_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_curve.md)
  or
  [`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md):
  one row per swept value, with the swept value, a decision rate, and
  the number of replicates behind it.

- target:

  The decision rate to solve for, strictly between 0 and 1.

- x:

  Name of the column holding the swept value. `NULL`, the default, takes
  the leading column.

- y:

  Name of the column holding the decision rate. `NULL`, the default,
  takes `power` or `p_meaningful`, whichever is present.

- n:

  The number of replicates behind each rate, either the name of a column
  or a numeric value. `NULL`, the default, takes `n_returned`,
  `n_converged` or `n_sims`, whichever is present.

- effect:

  Which focal effect to solve for, when the curve holds more than one.
  Matched against the `effect` or `param` column. `NULL`, the default,
  uses every row, which is correct only when the curve holds one effect.

- transform:

  The scale the swept value is fitted on: `"sqrt"` (the default, for a
  sample size), `"identity"` or `"log"`.

- level:

  Confidence level for the reported interval.

## Value

A list with elements `value` (the solved swept value), `lo` and `hi`
(its confidence bounds), `level`, `target`, `se` (the delta-method
standard error on the fitted scale, the scale on which the interval is
symmetric), `dispersion` (the heterogeneity factor applied, 1 where the
model fits), `x` and `y` (the columns used), `transform`, `intercept`
and `slope` (the fitted coefficients), `n_points` (the number of curve
points the fit used), and `x_min` and `x_max` (the swept range). A bound
is allowed to fall outside that range. When one does, the sweep was too
narrow to pin the value down and should be widened. A `dispersion` well
above 1 says the curve is not the shape the model assumes, so the solve
deserves a wider grid or more replicates before it deserves any trust.

## Details

The input is the data frame
[`power_curve_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_curve_mixed.md),
[`precision_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_curve.md)
or
[`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md)
returns, used as it stands. The swept value is taken from the leading
column, which is where all three put it, and the rate from `power` or
`p_meaningful`, whichever the curve carries. Each rate is a proportion
over a known number of replicates, and that count, read from
`n_returned`, `n_converged` or `n_sims`, weights the fit: a rate over
200 replicates should count for more than a rate over 20.

The fit is a binomial regression with a probit link, and the solved
value is the swept value at which the fitted rate equals `target`. The
probit is chosen because power is a normal tail probability: under the
normal approximation to a two-group comparison, the probit of power is
linear in the square root of the sample size, so the model has the shape
a design analysis already implies. Measured against
[`stats::power.t.test()`](https://rdrr.io/r/stats/power.t.test.html)
across twelve combinations of effect size and target power on each of
three grid shapes, at 400 replicates a point, the solved sample size
fell within 2.9% of the analytic answer on average, against 3.4% for a
logit fitted the same way.

The interval is the delta-method interval of
[`MASS::dose.p()`](https://rdrr.io/pkg/MASS/man/dose.p.html), computed
on the scale named by `transform` and mapped back, so it is symmetric on
that scale and asymmetric on the natural one. That asymmetry is the
honest shape: at the top of a power curve a given change in rate costs
far more sample size than the same change lower down. Where the
two-parameter model does not describe the curve, the interval is widened
by the heterogeneity factor of probit analysis, Pearson's chi-square
over its degrees of freedom (Finney, 1971), reported as `dispersion`. It
is floored at 1, so a well-fitting curve is left alone and a
badly-fitting one cannot report a narrower interval than its own
residuals justify.

The default `transform` of `"sqrt"` suits a sample-size axis, where a
rate rises with the square root of the sample size. Sweep something
else, an effect size or a random-effect standard deviation, and
`"identity"` is usually right.

Nothing here extrapolates. A curve whose rates do not straddle the
target is refused, with the range it did cover reported, and so is a fit
that solves outside the swept range. A curve whose fitted slope cannot
be told from zero is refused too: the crossing is then compatible with
any value at all, and an interval that said otherwise would be false.

What the interval covers is the Monte Carlo uncertainty of the fit, not
the gap between the fitted shape and the true curve. Across 36 checks
against
[`stats::power.t.test()`](https://rdrr.io/r/stats/power.t.test.html) at
400 replicates a point, the solved size sat within 2.9% of the analytic
answer on average and within 6.9% at worst, and the nominal 95% interval
covered the analytic value 35 times out of 36. Nearly all of that error
is the Monte Carlo noise the interval is describing, and it falls with
the square root of the replicate count. Raise the count far enough and
the interval narrows onto a fitted shape that is still slightly the
wrong shape, so replicates alone do not make a solved size arbitrarily
accurate. The remedies are a finer grid, more replicates, or a design
with a closed form to check against.

## References

Fieller, E. C. (1954). Some problems in interval estimation. *Journal of
the Royal Statistical Society: Series B*, 16(2), 175-185.
[doi:10.1111/j.2517-6161.1954.tb00159.x](https://doi.org/10.1111/j.2517-6161.1954.tb00159.x)

Finney, D. J. (1971). *Probit analysis* (3rd ed.). Cambridge University
Press.

## See also

[`target_n()`](https://pablobernabeu.github.io/pilotr/r/reference/target_n.md)
for the sample-size case, and
[`power_curve_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_curve_mixed.md),
[`precision_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_curve.md)
and
[`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md)
for the curves this consumes.

## Examples

``` r
spec <- build_spec(list(name = "s", seed = 1, design_kind = "between", n_subject = 40,
  factor_name = "group", lev1 = "a", lev2 = "b", intercept = 0, effect = 0.7,
  family = "gaussian", resp_name = "score", sigma = 1))
# n_sims is small so the example runs quickly. Use 200 or more for real planning.
curve <- sweep_spec(spec, "units$subject$n", c(20, 40, 60, 80), power_design, n_sims = 50)
solved <- solve_curve(curve, target = 0.8)
unlist(solved[c("value", "lo", "hi")])
#>    value       lo       hi 
#> 62.40150 52.15907 73.56151 

# This design has an analytic answer to check against, in total subjects across the two
# groups. It falls inside the interval, which fifty replicates a point make a wide one.
2 * stats::power.t.test(delta = 0.7, sd = 1, power = 0.8)$n
#> [1] 66.04934
```
