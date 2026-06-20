---
title: "Precision and ROPE design analysis"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Precision and ROPE design analysis}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




``` r
library(pilotr)
```

When sample sizes are large, or when the real question is whether an effect is *practically*
meaningful rather than merely non-zero, power against a point null is the wrong target. pilotr
implements a **precision-based design analysis** against a **region of practical equivalence**
(ROPE) — a fast, frequentist analog of the Bayesian highest-density-interval-versus-ROPE
approach.

> As in the power vignette, `n_sims` is kept tiny here so the vignette builds quickly; use
> `n_sims >= 200` for real planning.

## The idea

Over Monte Carlo replicates, pilotr records — for each focal fixed effect — the 95% confidence
interval and whether it falls **determinately outside** the ROPE (the effect is practically
meaningful) or **entirely inside** it (practical equivalence to zero), together with the
expected interval width. Sweeping sample size then locates the minimum N at which a focal
effect reaches a determinate decision with a target probability.

## A worked design

We reuse the crossed priming design — a small priming effect on log reaction time:


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

So that the analysis can run from the specification alone, pilotr auto-derives the analysis
model — the contrast columns, the response transform, and the mixed-model formula:


``` r
model_formula(spec_c)
#> .y ~ effect + (1 + effect | subject) + (1 + effect | item)
#> <environment: 0x0000021d71aa0f20>
```

## Precision at a fixed N

We declare the focal effect (its coefficient name and true value) and a ROPE half-width. Here
the effect lives on the log scale; we treat anything smaller than 0.02 as practically
equivalent to zero:


``` r
pr <- precision_design(spec_c, focal = c(effect = 0.05), rope = 0.02, n_sims = 25)
pr
#>    param true mean_ci_width p_meaningful p_equivalent n_converged
#> 1 effect 0.05    0.09491068          0.2            0          25
```

The columns: `p_meaningful` is the probability the 95% CI lands entirely outside the ROPE (a
determinate "meaningful" decision); `p_equivalent` is the probability it lands entirely inside
(a determinate "equivalent" decision); `mean_ci_width` is the expected precision.

## Sweeping sample size


``` r
prc <- precision_curve(spec_c, focal = c(effect = 0.05), subject_ns = c(15, 30), n_sims = 15)
prc[, c("n_subject", "param", "p_meaningful", "p_equivalent", "mean_ci_width")]
#>   n_subject  param p_meaningful p_equivalent mean_ci_width
#> 1        15 effect   0.00000000            0    0.12283252
#> 2        30 effect   0.06666667            0    0.08389208
```

As N grows the interval tightens, `p_meaningful` rises, and the design reaches a determinate
decision more reliably. Use this to choose the smallest N that meets a target — for example,
"`p_meaningful` ≥ 0.90" — rather than chasing power against a point null.

## See also

The [power vignette](power-analysis.html) covers simulation-based power and Type S/M errors;
the [getting-started vignette](getting-started.html) covers the core simulate-and-inspect loop.
