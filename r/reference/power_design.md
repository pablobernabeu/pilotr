# Simulation-based power and design analysis for a two-group Gaussian design

Estimate power by repeatedly simulating from the specification and
applying a two-sample t-test, alongside the Type S (sign) and Type M
(magnitude) design-analysis errors of Gelman and Carlin (2014).

## Usage

``` r
power_design(spec, n_sims = 1000, alpha = 0.05, workers = 1)
```

## Arguments

- spec:

  A design specification (path or list) for a two-group Gaussian design
  with one row per subject (see Details).

- n_sims:

  Number of Monte Carlo replicates. A power estimate carries a Monte
  Carlo standard error of about `sqrt(p * (1 - p) / n_sims)`, and
  `type_s` and `type_m` average over the significant replicates alone,
  so they settle more slowly still. At least 200 replicates are
  advisable for study planning.

- alpha:

  Two-sided significance level.

- workers:

  Number of local worker processes over which to spread the replicates.
  The default of 1 runs serially. Because every replicate takes its own
  seed from
  [`replicate_seeds()`](https://pablobernabeu.github.io/pilotr/r/reference/replicate_seeds.md),
  any worker count returns results identical to a serial run.

## Value

A list with elements `n_sims`, `alpha`, `power`, `n_significant`,
`true_effect`, `mean_estimate`, `type_s` (sign-error rate among
significant replicates), and `type_m` (mean exaggeration ratio among
significant replicates). Both design-analysis quantities are `NaN` when
no replicate reached significance and when the true effect is zero, as
in the null condition
[`design_conditions()`](https://pablobernabeu.github.io/pilotr/r/reference/design_conditions.md)
produces: neither is defined without a true value to compare against,
and Type M divides by it.

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

A specification with no coefficient for the factor's contrast has a true
effect of 0, so `type_s` and `type_m` are `NaN`.

## References

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
#> $power
#> [1] 0.52
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
#> $type_m
#> [1] 1.303814
#> 
```
