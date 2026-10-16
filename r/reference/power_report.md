# Write a sample-size justification from a power result

Turns the result of
[`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md),
[`power_design()`](https://pablobernabeu.github.io/pilotr/r/reference/power_design.md)
or
[`precision_design()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md),
with the
[`target_n()`](https://pablobernabeu.github.io/pilotr/r/reference/target_n.md)
result that solved a curve for the sample size, into a sample-size
justification in Markdown, ready for a preregistration, a registered
report or a grant application.

## Usage

``` r
power_report(
  x,
  solved = NULL,
  spec_file = NULL,
  justification = NULL,
  null_result = NULL,
  rope = NULL
)

# S3 method for class 'pilotr_report'
print(x, ...)
```

## Arguments

- x:

  For `power_report()`, the result of
  [`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md),
  [`power_design()`](https://pablobernabeu.github.io/pilotr/r/reference/power_design.md)
  or
  [`precision_design()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md)
  at one design. A curve goes to
  [`target_n()`](https://pablobernabeu.github.io/pilotr/r/reference/target_n.md),
  and its result to `solved`. For
  [`print()`](https://rdrr.io/r/base/print.html), a report.

- solved:

  Optional result of
  [`target_n()`](https://pablobernabeu.github.io/pilotr/r/reference/target_n.md),
  or of
  [`solve_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/solve_curve.md)
  for an axis other than sample size, which gives the desired power and
  the value solved for it.

- spec_file:

  Optional path of the specification file the analysis read. The report
  gives the MD5 checksum of its bytes, its seed and its units.

- justification:

  Optional justification of each focal effect's size: one string for
  every focal effect, or a character vector named by focal effect. An
  effect without one is reported as having none.

- null_result:

  Optional result of the same analysis with every focal effect set to
  zero. A result whose focal effects are not zero is refused.

- rope:

  The region of practical equivalence given to
  [`precision_design()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md),
  as a half-width or an interval, which the report states. It applies to
  a precision result alone.

- ...:

  Ignored, present for consistency with the generic.

## Value

A character string of class `pilotr_report` holding the Markdown, which
ends in a newline. Printing it shows the text, and
`cat(report, file = "justification.md")` saves it. The same arguments
always give the same text.

## Details

The sections follow the recommendations of Lakens (2022, Table 4) for
reporting an a-priori power analysis. The primary analysis comes first,
with its model, test, alpha and sidedness, then the desired power and
the sample size solved for it, with its interval. Each focal effect
follows, with its true value, the justification given for that value and
its power. Last comes the power when the null hypothesis is true. The
simulation is then accounted for as Morris, White and Crowther (2019)
recommend for simulation studies: the replicates simulated and those
that returned an estimate, the fits that were singular or carried a
warning, and the Monte Carlo standard error of every rate. The report
ends with what a reader needs to reproduce the run: the software and its
versions, the specification's seed and an MD5 checksum of the
specification file's bytes. The versions are those loaded when the
report is written, so write it in the session that ran the analysis.

Each rate's Monte Carlo standard error and 95% Wilson interval are
computed in the report from the number of significant replicates and the
number of replicates behind the rate, with the functions that compute
them in the results. Type S has its Wilson interval from the number of
wrong signs among the significant replicates. The Python twin's
`power_report()` writes the same text, byte for byte, from the same
counts, which `tools/parity/report_cross.py` in the source repository
checks.

The rates are those of the design `x` simulated. A sample size solved
from a curve need not be that design's, so a justification usually pairs
the solved size with a run of
[`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md)
at that size. `spec_file` names the specification of that run, and the
report lists its units. A result carries neither its specification nor
its seed, so without `spec_file` the report says that both are missing.

`null_result` is the result of the same analysis of a design whose focal
effects are all zero, such as the null condition that
[`design_conditions()`](https://pablobernabeu.github.io/pilotr/r/reference/design_conditions.md)
puts first. Its rejection rate is the Type I error rate of the test,
which for a Wald z with few subjects or items can exceed alpha, and
which a preregistered simplification rule can raise (see
[`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md)).

A
[`precision_design()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md)
result gives a report on the interval and the region of practical
equivalence in place of the test. Each focal effect then has its mean
interval width, the probability of a meaningful effect, with the
interval wholly outside the region, and the probability of practical
equivalence, with the interval wholly inside it. `rope` states the
region, since the result does not record it. A
[`target_n()`](https://pablobernabeu.github.io/pilotr/r/reference/target_n.md)
result solved over `p_meaningful` gives the desired probability of a
meaningful decision. This section exists in R alone, as
[`precision_design()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md)
does.

## References

Lakens, D. (2022). Sample size justification. *Collabra: Psychology*,
8(1), 33267.
[doi:10.1525/collabra.33267](https://doi.org/10.1525/collabra.33267)

Morris, T. P., White, I. R. and Crowther, M. J. (2019). Using simulation
studies to evaluate statistical methods. *Statistics in Medicine*,
38(11), 2074-2102.
[doi:10.1002/sim.8086](https://doi.org/10.1002/sim.8086)

## See also

[`target_n()`](https://pablobernabeu.github.io/pilotr/r/reference/target_n.md)
for the sample size a curve implies, and
[`design_conditions()`](https://pablobernabeu.github.io/pilotr/r/reference/design_conditions.md)
for the null condition.

## Examples

``` r
spec <- build_spec(list(name = "d", seed = 1, design_kind = "between",
  factor_name = "group", lev1 = "a", lev2 = "b", n_subject = 64,
  intercept = 100, effect = 5, family = "gaussian", resp_name = "", sigma = 10))
# n_sims is small so the example runs quickly. Use 200 or more for real planning.
res <- power_design(spec, n_sims = 50)
spec_file <- tempfile(fileext = ".json")
writeLines(spec_json(spec), spec_file)
power_report(res, spec_file = spec_file,
  justification = "Half a standard deviation, the smallest difference worth acting on.")
#> # Sample-size justification
#> 
#> ## Primary analysis
#> 
#> - Model and test: stats::t.test, two-sample t with pooled variance
#> - Significance level: alpha = 0.05, two-sided
#> - Simulated units: subject 64
#> 
#> ## Focal effects
#> 
#> The rates below are those of the simulated design, each over the replicates that returned an estimate for the effect.
#> 
#> ### The difference between the two groups
#> 
#> - True value: 5
#> - Justification: Half a standard deviation, the smallest difference worth acting on.
#> - Power: 0.520 (Monte Carlo standard error 0.071; 95% Wilson interval 0.385 to 0.652), from 26 of 50 replicates
#> - Type S: 0.000 (95% Wilson interval 0.000 to 0.129)
#> - Type M: 1.304 (Monte Carlo standard error 0.047)
#> 
#> ## Monte Carlo accounting
#> 
#> - Replicates: 50 simulated
#> 
#> Each rate's Monte Carlo standard error is sqrt(p(1 - p) / n), over the n replicates that returned an estimate, and its interval is the 95% Wilson score interval. The standard error of Type M is the standard deviation of the exaggeration ratios over the square root of their number.
#> 
#> ## Reproducibility
#> 
#> - Software: pilotr 0.3.1, R 4.6.1
#> - Specification: file19fa75d164f8.json, MD5 checksum fb7cf3e60def2bacdd2e7508dc1575be
#> - Seed: 1
```
