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

  A design specification (path or list) for a two-group Gaussian design.

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
  The default of 1 runs serially. Because every replicate seeds the
  shared RNG from its own index, any worker count returns results
  identical to a serial run.

## Value

A list with elements `n_sims`, `alpha`, `power`, `n_significant`,
`true_effect`, `mean_estimate`, `type_s` (sign-error rate among
significant replicates), and `type_m` (mean exaggeration ratio among
significant replicates).

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
#> [1] 0.6
#> 
#> $n_significant
#> [1] 30
#> 
#> $true_effect
#> [1] 5
#> 
#> $mean_estimate
#> [1] 5.238285
#> 
#> $type_s
#> [1] 0
#> 
#> $type_m
#> [1] 1.361989
#> 
```
