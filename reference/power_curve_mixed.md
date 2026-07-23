# Power curve over sample size for a crossed mixed-effects design

Sweep the number of subjects and compute mixed-effects power at each, to
locate where power crosses a target. Runs the same replicate loop as
[`power_mixed()`](https://pablobernabeu.github.io/pilotr/reference/power_mixed.md)
and so requires the `lme4` and `lmerTest` packages.

## Usage

``` r
power_curve_mixed(spec, subject_ns, n_sims = 60, alpha = 0.05, workers = 1)
```

## Arguments

- spec:

  A design specification (path or list) with one within-unit factor.

- subject_ns:

  A numeric vector of subject counts to evaluate.

- n_sims:

  Number of Monte Carlo replicates per point. A power estimate carries a
  Monte Carlo standard error of about `sqrt(p * (1 - p) / n_sims)`, and
  `type_m` averages over the significant replicates alone, so it settles
  more slowly still. At least 200 replicates are advisable for study
  planning.

- alpha:

  Two-sided significance level.

- workers:

  Number of local worker processes over which to spread the replicates
  at each grid point. The default of 1 runs serially, and any worker
  count returns results identical to a serial run.

## Value

A data frame with one row per sample size and columns `n_subject`,
`power`, `type_m`, and `n_converged`. As in
[`power_mixed()`](https://pablobernabeu.github.io/pilotr/reference/power_mixed.md),
`power` at each sample size is the proportion of significant results
among the `n_converged` converged replicates.

## Details

Like
[`power_mixed()`](https://pablobernabeu.github.io/pilotr/reference/power_mixed.md),
this function runs pilotr's own simulation loop over the portable design
specification rather than wrapping an existing package. It covers
territory pioneered by `simr` (Green and MacLeod, 2016) and `mixedpower`
(Kumle, Vo and Draschkow, 2021), and differs in being driven by the
portable cross-language specification, in reporting Type S and Type M
errors, and in built-in parallelisation: with `workers > 1` a single
worker pool is created once and reused across all sample sizes in the
sweep.

## References

Green, P. and MacLeod, C. J. (2016). SIMR: An R package for power
analysis of generalized linear mixed models by simulation. *Methods in
Ecology and Evolution*, 7(4), 493-498.
[doi:10.1111/2041-210x.12504](https://doi.org/10.1111/2041-210x.12504)

Kumle, L., Vo, M. L.-H. and Draschkow, D. (2021). Estimating power in
(generalized) linear mixed models: An open introduction and tutorial in
R. *Behavior Research Methods*, 53, 2528-2543.
[doi:10.3758/s13428-021-01546-0](https://doi.org/10.3758/s13428-021-01546-0)

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
  # n_sims is small so the example runs quickly. Use 200 or more for real planning.
  power_curve_mixed(spec, subject_ns = c(12, 18), n_sims = 8)
}
#>   n_subject power   type_m n_converged
#> 1        12  0.25 2.573080           8
#> 2        18  0.25 2.489337           8
# }
```
