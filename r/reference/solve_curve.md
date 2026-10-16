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
  the leading column, which holds a grid index when the sweep was over
  lists.

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
fell within 2.6% of the analytic answer on average. A logit fitted the
same way did as well on that draw, at 2.5%, and worse on another, at
3.5% against 2.2%. Curves this noisy do not separate the two links, and
the probit is chosen for its shape.

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

A sweep over the conditions
[`design_conditions()`](https://pablobernabeu.github.io/pilotr/r/reference/design_conditions.md)
builds is a sweep over lists, so its leading column holds a grid index,
and solving it as it stands would solve for that index. Either name the
column holding the effect with `x`, such as `true` in the rows of
[`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md),
or sweep one coefficient by its own path. The curve from
`sweep_spec(spec, "fixed$coefficients$cond", c(0, 0.02, 0.04), power_mixed)`
holds the other coefficients at their values and leads with a column
named `cond`. From that curve,
`solve_curve(curve, target = 0.8, effect = "cond", transform = "identity")`
returns the minimum detectable effect.

Nothing here extrapolates. A curve whose rates do not straddle the
target is refused, with the range it did cover reported, and so is a fit
that solves outside the swept range. A curve whose fitted slope cannot
be told from zero is refused too: the crossing is then compatible with
any value at all, and an interval that said otherwise would be false.

A curve over subjects in a crossed design need not rise towards 1. When
an effect varies from item to item, as a by-item random slope says it
does, power at a fixed number of items levels off below 1 (Westfall,
Kenny and Judd, 2014; Judd, Westfall and Kenny, 2017). Adding subjects
cannot lift it past that ceiling. The probit model rises towards 1 and
cannot follow such a plateau. A curve that levels off near the target
can then be solved far from the true crossing. One that levels off just
below the target never crosses it, but Monte Carlo noise can lift a
simulated rate over the target, and the curve can then be solved as
though it crossed. To compare designs that trade subjects for items,
solve curves over the number of items, one at each of a few subject
counts, from
[`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md)
over `units$item$n`. Power over items levels off in the same way, at a
ceiling set by the by-subject slope and the number of subjects. The
subject counts must be large enough to put that ceiling well above the
target.

What the interval covers is the Monte Carlo uncertainty of the fit, not
the gap between the fitted shape and the true curve, and it treats the
points of the curve as independent. Every sweep in pilotr now draws each
point's replicates from a seed of its own (see
[`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md)).
Over 1,000 simulated curves for each combination of four effect sizes
and three target powers, at 400 replicates a point, the nominal 95%
interval covered the analytic sample size in 96.5%, 95.9% and 96.2% of
solves on the three grid shapes. For an effect size swept at a fixed
sample size, it covered the analytic effect in 95.7% of solves. When
every point reused the specification's seed, as up to 0.3.1, the figures
were 87.8%, 94.0%, 93.5% and 70.3%. A curve assembled by hand, from
separate runs that each started from the same seed, shares its noise in
the same way. In the calibration's 36 checks, the solved size sat within
2.6% of the analytic answer on average and within 7.7% at worst. Nearly
all of that error is the Monte Carlo noise the interval is describing,
and it falls with the square root of the replicate count. Raise the
count far enough and the interval narrows onto a fitted shape that is
still slightly the wrong shape, so replicates alone do not make a solved
size arbitrarily accurate. The remedies are more replicates, a finer
grid or a design with a closed form to check against. A finer grid helps
only while its points are independent. With shared seeds, its points
share more of their noise, and the tight grid, whose points lie closest
together, covered 87.8% where the coarse one covered 94.0%.

## References

Fieller, E. C. (1954). Some problems in interval estimation. *Journal of
the Royal Statistical Society: Series B*, 16(2), 175-185.
[doi:10.1111/j.2517-6161.1954.tb00159.x](https://doi.org/10.1111/j.2517-6161.1954.tb00159.x)

Finney, D. J. (1971). *Probit analysis* (3rd ed.). Cambridge University
Press.

Judd, C. M., Westfall, J. and Kenny, D. A. (2017). Experiments with more
than one random factor: Designs, analytic models, and statistical power.
*Annual Review of Psychology*, 68, 601-625.
[doi:10.1146/annurev-psych-122414-033702](https://doi.org/10.1146/annurev-psych-122414-033702)

Westfall, J., Kenny, D. A. and Judd, C. M. (2014). Statistical power and
optimal design in experiments in which samples of participants respond
to samples of stimuli. *Journal of Experimental Psychology: General*,
143(5), 2020-2045.
[doi:10.1037/xge0000014](https://doi.org/10.1037/xge0000014)

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
#> 70.48928 56.16037 86.44456 

# This design has an analytic answer to check against, in total subjects across the two
# groups. It falls inside the interval, which fifty replicates a point make a wide one.
2 * stats::power.t.test(delta = 0.7, sd = 1, power = 0.8)$n
#> [1] 66.04934
```
