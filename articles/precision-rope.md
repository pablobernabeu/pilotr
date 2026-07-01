# Precision and ROPE design analysis

``` r

library(pilotr)
```

When sample sizes are large, or when the question of interest is whether
an effect is *practically* meaningful, power against a point null is
often the wrong target. pilotr implements a precision-based design
analysis against a region of practical equivalence (ROPE). This is a
fast, frequentist analogue of the Bayesian approach that compares a
highest-density interval with a ROPE.

> As in the power vignette, `n_sims` is kept small here so the vignette
> builds quickly. For real planning, a value of `n_sims >= 200` is
> recommended.

## The idea

Over Monte Carlo replicates, and for each focal fixed effect, pilotr
records the 95% confidence interval and whether it falls determinately
outside the ROPE (the effect is practically meaningful) or entirely
inside it (practical equivalence to zero), together with the expected
interval width. Sweeping sample size then locates the minimum *N* at
which a focal effect reaches a determinate decision with a target
probability.

## A worked design

We reuse the crossed priming design, which represents a small priming
effect on log reaction time.

``` r

spec_c <- build_spec(list(
  name = "priming", seed = 1, design_kind = "within", include_items = TRUE,
  n_subject = 24, n_item = 18,
  factor_name = "condition", lev1 = "related", lev2 = "unrelated",
  intercept = 6, effect = 0.05,
  subj_int_sd = 0.12, subj_slope_sd = 0.04, subj_corr = 0.2,
  item_int_sd = 0.08, item_slope_sd = 0.02, item_corr = -0.1,
  family = "shifted_lognormal", resp_name = "RT", sigma = 0.3, shift = 200))
```

So that the analysis can run from the specification alone, pilotr
auto-derives the analysis model, including the contrast columns, the
response transform and the mixed-model formula.

``` r

model_formula(spec_c)
#> .y ~ effect + (1 + effect | subject) + (1 + effect | item)
#> <environment: 0x55b3252b2c38>
```

## Precision at a fixed *N*

We declare the focal effect (its coefficient name and true value) and a
ROPE half-width. Here the effect lives on the log scale, and we treat
anything smaller than 0.02 as practically equivalent to zero.

``` r

pr <- precision_design(spec_c, focal = c(effect = 0.05), rope = 0.02, n_sims = 25)
pr
#>    param true mean_ci_width p_meaningful p_equivalent n_converged
#> 1 effect 0.05    0.09491084          0.2            0          25
```

The columns are interpreted as follows. `p_meaningful` is the
probability the 95% CI lands entirely outside the ROPE, a determinate
‘meaningful’ decision. `p_equivalent` is the probability it lands
entirely inside, a determinate ‘equivalent’ decision. `mean_ci_width` is
the expected precision.

## Sweeping sample size

``` r

prc <- precision_curve(spec_c, focal = c(effect = 0.05), subject_ns = c(15, 30), n_sims = 15)
prc[, c("n_subject", "param", "p_meaningful", "p_equivalent", "mean_ci_width")]
#>   n_subject  param p_meaningful p_equivalent mean_ci_width
#> 1        15 effect   0.00000000            0    0.12283251
#> 2        30 effect   0.06666667            0    0.08389208
```

As *N* grows the interval tightens, `p_meaningful` rises and the design
reaches a determinate decision more reliably. This curve can be used to
choose the smallest *N* that meets a target, for example ‘`p_meaningful`
≥ 0.90’, which provides a more informative criterion than power against
a point null.

## See also

The [power
vignette](https://pablobernabeu.github.io/pilotr/articles/power-analysis.md)
covers simulation-based power and Type S/M errors. The [getting-started
vignette](https://pablobernabeu.github.io/pilotr/articles/getting-started.md)
covers the core simulate-and-inspect loop.
