# Precision and ROPE design analysis at a fixed sample size

A fast frequentist analogue of a Bayesian
highest-density-interval-versus-ROPE design analysis. Across Monte Carlo
replicates, fit the model and record, for each focal fixed effect,
whether its 95% confidence interval falls entirely outside a region of
practical equivalence (a practically meaningful effect) or entirely
inside it (practical equivalence to zero), along with the expected
interval width. Requires the `lme4` package.

## Usage

``` r
precision_design(
  spec,
  focal = NULL,
  formula = NULL,
  prep = NULL,
  rope = 0.05,
  n_sims = 100,
  workers = 1
)
```

## Arguments

- spec:

  A design specification (path or list).

- focal:

  The focal effects. `NULL`, the default, analyses every coefficient in
  the specification and takes the true values from it. A named numeric
  vector maps coefficient names to their true values, and a character
  vector names them without their true values. Interaction effects
  follow the model's column naming, so a specification key `a:b` is the
  focal name `a_b`.

- formula:

  Optional `lme4` formula; if `NULL` it is derived from the
  specification via
  [`model_formula()`](https://pablobernabeu.github.io/pilotr/r/reference/model_formula.md).

- prep:

  Optional function mapping a simulated data set to the modelling data;
  if `NULL` it is derived via
  [`model_data()`](https://pablobernabeu.github.io/pilotr/r/reference/model_data.md),
  which log-transforms the outcome and builds the contrast and
  interaction columns, so focal names follow the auto-formula
  (interactions written as `a_b`).

- rope:

  Half-width of the region of practical equivalence; an effect with
  `abs(beta) < rope` is treated as practically equivalent to zero. Set
  it clearly narrower than the smallest effect worth detecting, because
  the probability of a determinate meaningful decision about an effect
  no larger than `rope` cannot rise above 0.5 however large the sample.

- n_sims:

  Number of Monte Carlo replicates. `p_meaningful` and `p_equivalent`
  are proportions over the converged replicates, so they carry a Monte
  Carlo standard error of about `sqrt(p * (1 - p) / n_sims)` and move in
  coarse steps when `n_sims` is small. At least 200 replicates are
  advisable for real planning.

- workers:

  Number of local worker processes over which to spread the replicates.
  The default of 1 runs serially. Because every replicate seeds the
  shared RNG from its own index, any worker count returns results
  identical to a serial run.

## Value

A data frame with one row per focal effect and columns `param`, `true`,
`mean_ci_width`, `p_meaningful`, `p_equivalent`, `n_attempted`,
`n_returned`, `n_converged`, `n_singular`, and `n_warning`. The interval
behind `mean_ci_width` and the ROPE decisions is the Wald approximation
described in Details.

The decision proportions are taken over `n_returned`, the replicates
that produced an estimate. The remaining counts separate the fit
outcomes, because a fit can return a usable estimate while still being
boundary-singular or carrying a convergence warning: `n_converged`
counts replicates with neither, `n_singular` those where
[`lme4::isSingular()`](https://rdrr.io/pkg/lme4/man/isSingular.html) was
true, and `n_warning` those with a warning or optimiser convergence
message. Singular and warning fits are retained, since their
fixed-effect estimates remain interpretable and discarding them would
bias the result. A large `n_singular` means the model being fitted is
richer than the design can support at that sample size, which is common
in crossed designs (Bates et al., 2015; Matuschek et al., 2017).

## Details

The interval is a Wald approximation: the estimate plus or minus 1.96
standard errors from the model's variance-covariance matrix. This
fixed-z interval is chosen for speed and for comparability across
replicates; in small samples it is somewhat narrower than a
Satterthwaite t interval, so `p_meaningful` and `mean_ci_width` are
slightly optimistic at small sample sizes.

## References

Bates, D., Kliegl, R., Vasishth, S. and Baayen, H. (2015). Parsimonious
mixed models. *arXiv*.
[doi:10.48550/arXiv.1506.04967](https://doi.org/10.48550/arXiv.1506.04967)

Matuschek, H., Kliegl, R., Vasishth, S., Baayen, H. and Bates, D.
(2017). Balancing Type I error and power in linear mixed models.
*Journal of Memory and Language*, 94, 305-315.
[doi:10.1016/j.jml.2017.01.001](https://doi.org/10.1016/j.jml.2017.01.001)

## Examples

``` r
# \donttest{
if (requireNamespace("lme4", quietly = TRUE)) {
  spec <- build_spec(list(name = "pr", seed = 1, design_kind = "within",
    include_items = TRUE, n_subject = 12, n_item = 12, factor_name = "cond",
    lev1 = "a", lev2 = "b", intercept = 6, effect = 0.05,
    subj_int_sd = 0.12, subj_slope_sd = 0.04, subj_corr = 0.2,
    item_int_sd = 0.08, item_slope_sd = 0.02, item_corr = -0.1,
    family = "shifted_lognormal", resp_name = "", sigma = 0.3, shift = 200))
  # n_sims is small so the example runs quickly. Use 200 or more for real planning.
  precision_design(spec, focal = c(effect = 0.05), rope = 0.02, n_sims = 10)
}
#>    param true mean_ci_width p_meaningful p_meaningful_mcse p_meaningful_lo
#> 1 effect 0.05     0.1526052          0.1        0.09486833      0.01787621
#>   p_meaningful_hi p_equivalent p_equivalent_mcse p_equivalent_lo
#> 1         0.40415            0                 0               0
#>   p_equivalent_hi n_attempted n_returned n_converged n_singular n_warning
#> 1       0.2775328          10         10           0         10        10
# }
```
