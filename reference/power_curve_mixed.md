# Power curve over sample size for a crossed mixed-effects design

Sweep the number of subjects and compute mixed-effects power at each, to
locate where power crosses a target. Calls \[power_mixed()\] and so
requires the \`lme4\` and \`lmerTest\` packages.

## Usage

``` r
power_curve_mixed(spec, subject_ns, n_sims = 60, alpha = 0.05)
```

## Arguments

- spec:

  A design specification (path or list) with one within-unit factor.

- subject_ns:

  A numeric vector of subject counts to evaluate.

- n_sims:

  Number of Monte Carlo replicates per point.

- alpha:

  Two-sided significance level.

## Value

A data frame with one row per sample size and columns \`n_subject\`,
\`power\`, and \`type_m\`.

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
  power_curve_mixed(spec, subject_ns = c(12, 18), n_sims = 8)
}
#>   n_subject power   type_m
#> 1        12  0.25 2.573080
#> 2        18  0.25 2.489337
# }
```
