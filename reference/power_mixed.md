# Simulation-based power and design analysis for a crossed mixed-effects design

For each replicate, simulate from the ground-truth specification, fit
the maximal model \`y ~ cond + (1 + cond \| subject) + (1 + cond \|
item)\` with \`lmerTest\`, and test the fixed effect using Satterthwaite
p-values. Reports power together with the Type S and Type M errors of
Gelman and Carlin (2014). Requires the \`lme4\` and \`lmerTest\`
packages.

## Usage

``` r
power_mixed(spec, n_sims = 100, alpha = 0.05)
```

## Arguments

- spec:

  A design specification (path or list) with exactly one within-unit
  factor.

- n_sims:

  Number of Monte Carlo replicates.

- alpha:

  Two-sided significance level.

## Value

A list with elements \`n_sims\`, \`n_converged\`, \`alpha\`, \`power\`,
\`n_significant\`, \`true_effect\`, \`mean_estimate\`, \`type_s\`, and
\`type_m\`.

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
  power_mixed(spec, n_sims = 10)
}
#> $n_sims
#> [1] 10
#> 
#> $n_converged
#> [1] 10
#> 
#> $alpha
#> [1] 0.05
#> 
#> $power
#> [1] 0.2
#> 
#> $n_significant
#> [1] 2
#> 
#> $true_effect
#> [1] 0.05
#> 
#> $mean_estimate
#> [1] 0.05286316
#> 
#> $type_s
#> [1] 0
#> 
#> $type_m
#> [1] 2.57308
#> 
# }
```
