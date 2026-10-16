# Simulation-based power and design analysis for a two-group Gaussian design

Estimate power by repeatedly simulating from the specification and
applying a two-sample t-test, alongside the Type S (sign) and Type M
(magnitude) design-analysis errors of Gelman and Carlin (2014).

## Usage

``` r
power_design(
  spec,
  n_sims = 1000,
  alpha = 0.05,
  workers = 1,
  cl = NULL,
  alternative = c("two.sided", "greater", "less")
)
```

## Arguments

- spec:

  A design specification (path or list) for a two-group Gaussian design
  with one row per subject (see Details).

- n_sims:

  Number of Monte Carlo replicates. A power estimate carries a Monte
  Carlo standard error of about `sqrt(p * (1 - p) / n_sims)`, reported
  as `power_mcse`, and `type_s` and `type_m` average over the
  significant replicates alone, so they settle more slowly still. At
  least 200 replicates are advisable for study planning.

- alpha:

  Significance level of the test, two-sided unless `alternative` says
  otherwise.

- workers:

  Number of local worker processes over which to spread the replicates.
  The default of 1 runs serially. Because every replicate takes its own
  seed from
  [`replicate_seeds()`](https://pablobernabeu.github.io/pilotr/r/reference/replicate_seeds.md),
  any worker count returns results identical to a serial run.

- cl:

  A cluster made by
  [`parallel::makeCluster()`](https://rdrr.io/r/parallel/makeCluster.html)
  to run the replicates on, in place of one started for this call, in
  which case `workers` is ignored. The cluster is left running for its
  owner to stop, and through
  [`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md)
  one cluster serves every grid point.

- alternative:

  The test's alternative hypothesis: `"two.sided"`, the default, or
  `"greater"` or `"less"` for a one-sided test of a directional
  prediction (see Details).

## Value

A list with elements `n_sims`, `alpha`, `alternative`, `power`,
`power_mcse` (its Monte Carlo standard error), `power_lo` and `power_hi`
(its 95% Wilson interval), `n_significant`, `true_effect` and
`mean_estimate`. Then come `type_s` (sign-error rate among significant
replicates), `type_s_lo` and `type_s_hi` (its 95% Wilson interval),
`type_m` (mean exaggeration ratio among significant replicates) and
`type_m_mcse` (its Monte Carlo standard error). Both design-analysis
quantities are `NaN` when no replicate reached significance and when the
true effect is zero, as in the null condition
[`design_conditions()`](https://pablobernabeu.github.io/pilotr/r/reference/design_conditions.md)
produces: neither is defined without a true value to compare against,
and Type M divides by it. Their interval and standard error are `NaN`
then too, and `type_m_mcse` is also `NaN` when a single replicate
reached significance, since a standard deviation needs two.

## Details

The t-test treats every row as an independent observation, which is
valid only when each subject contributes one row and no rows share a
cluster. `power_design()` therefore takes a Gaussian design with exactly
one two-level factor between subjects and no item unit, within factor or
grouping factor besides `subject`. Continuous predictors and by-subject
random effects are allowed, since with one row per subject they vary
independently from row to row. Any other design is refused with a
message naming what makes its rows correlated.
[`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md)
fits the model that such a design implies. Take 30 subjects crossed with
20 items, with by-subject and residual standard deviations of 1 and a
by-item one of 0.3. With no true effect, a t-test of every row is
significant in about half of all replicates.

Each level of the factor needs at least two subjects, so a design needs
at least 4. The subjects are split into two blocks in level order, and
the second level holds `floor(N / 2)`. With one subject at a level, the
pooled variance rests on the other level alone, and with one at each it
has no degree of freedom at all.

A specification with no coefficient for the factor's contrast has a true
effect of 0, so `type_s` and `type_m` are `NaN`.

Power, Type S and Type M are each estimated from a finite number of
replicates, so each is reported with its Monte Carlo error. Power has
its standard error and a 95% Wilson interval, one of the intervals
Brown, Cai and DasGupta (2001) recommend in place of the Wald interval.
It stays between 0 and 1 and does not shrink to a point when power is 0
or 1. Type S and Type M rest on the significant replicates alone, often
a handful. Type S has a Wilson interval, because its standard error is
zero whenever no sign is wrong. Type M has its standard error, the
standard deviation of the exaggeration ratios over the square root of
their number.

A directional prediction calls for a one-sided test, which `alternative`
selects. The estimate is the mean of the second level minus that of the
first, whose true value is `true_effect`. `"greater"` tests whether it
is positive and `"less"` whether it is negative, by
[`stats::t.test()`](https://rdrr.io/r/stats/t.test.html) with the same
`alternative`. The one-sided p-value is half the two-sided one when the
estimate lies in the tested direction. A one-sided test at `alpha`
therefore rejects in exactly the replicates in which the two-sided test
at `2 * alpha` rejects with the estimate in that direction. Type S and
Type M are taken over the replicates significant in the tested
direction. When that is the direction of the true effect, Type S is 0 by
construction, and Type M is the exaggeration among the replicates that
would support the prediction.

Results from earlier versions were two-sided. A one-sided power at
`alpha`, in the direction of the true effect, can be read off a
two-sided result at `2 * alpha`: it is `power * (1 - type_s)`. Type M is
recovered only when that result's Type S is 0, since the wrong-signed
replicates cannot otherwise be told apart from the others.

## References

Brown, L. D., Cai, T. T. and DasGupta, A. (2001). Interval estimation
for a binomial proportion. *Statistical Science*, 16(2), 101-133.
[doi:10.1214/ss/1009213286](https://doi.org/10.1214/ss/1009213286)

Gelman, A. and Carlin, J. (2014). Beyond power calculations: Assessing
Type S (sign) and Type M (magnitude) errors. *Perspectives on
Psychological Science*, 9(6), 641-651.
[doi:10.1177/1745691614551642](https://doi.org/10.1177/1745691614551642)

## Examples

``` r
spec <- build_spec(list(name = "d", seed = 1, design_kind = "between",
  factor_name = "group", lev1 = "a", lev2 = "b", n_subject = 64,
  intercept = 100, effect = 5, family = "gaussian", resp_name = "", sigma = 10))
# n_sims is small so the example runs quickly. Use 200 or more for real planning.
power_design(spec, n_sims = 50)
#> $n_sims
#> [1] 50
#> 
#> $alpha
#> [1] 0.05
#> 
#> $alternative
#> [1] "two.sided"
#> 
#> $power
#> [1] 0.52
#> 
#> $power_mcse
#> [1] 0.07065409
#> 
#> $power_lo
#> [1] 0.3851174
#> 
#> $power_hi
#> [1] 0.6520286
#> 
#> $n_significant
#> [1] 26
#> 
#> $true_effect
#> [1] 5
#> 
#> $mean_estimate
#> [1] 4.801842
#> 
#> $type_s
#> [1] 0
#> 
#> $type_s_lo
#> [1] 0
#> 
#> $type_s_hi
#> [1] 0.1287289
#> 
#> $type_m
#> [1] 1.303814
#> 
#> $type_m_mcse
#> [1] 0.04732007
#> 
```
