---
title: "Response families"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Response families}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




``` r
library(pilotr)
```

Behavioral outcomes are rarely Gaussian. pilotr maps the linear predictor
$\eta = \beta_0 + \sum_k \beta_k x_k + (\text{random effects})$ to an outcome through one of
six response families. Crucially, the **fixed intercept and effect live on the family's own
scale** — identity for Gaussian, the log scale for reaction times and counts, and the logit
scale for accuracy, ordinal, and proportion outcomes. This vignette shows each family with
scale-appropriate parameters.

A small helper to build a two-group between-subjects design and report the group means:


``` r
demo <- function(family, intercept, effect, n = 4000, ...) {
  spec <- build_spec(list(
    name = family, seed = 1, design_kind = "between", n_subject = n,
    factor_name = "group", lev1 = "control", lev2 = "treatment",
    intercept = intercept, effect = effect, family = family, resp_name = "", ...))
  d <- simulate_design(spec)
  y <- d[[spec$response$name]]
  list(spec = spec, data = d, y = y, by_group = tapply(y, d$group, mean))
}
```

## Gaussian

The default: $y = \eta + \varepsilon$, with residual standard deviation `sigma`. Use it for
continuous, roughly symmetric outcomes (ratings averaged over many trials, standardized
scores).


``` r
g <- demo("gaussian", intercept = 100, effect = 5, sigma = 10)
round(g$by_group, 2)
#>   control treatment 
#>     97.50    102.41
hist(g$y, breaks = 40, col = "#2C6FB0", border = "white", main = "Gaussian", xlab = "score")
```

![](C:/Users/pablob/AppData/Local/Temp/RtmpmW4XlM/response-families_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

## Shifted lognormal (reaction times)

Reaction times are right-skewed and bounded below. pilotr models them as
$y = \text{shift} + \exp(\eta + \varepsilon)$, with the effect on the **log** scale. An
intercept of 6 implies a typical RT near `exp(6)` ms above the shift.


``` r
rt <- demo("shifted_lognormal", intercept = 6, effect = 0.1, sigma = 0.3, shift = 200)
round(rt$by_group, 1)
#>   control treatment 
#>     600.6     642.6
hist(rt$y, breaks = 40, col = "#B0402C", border = "white", main = "Shifted lognormal (RT)",
     xlab = "RT (ms)")
```

![](C:/Users/pablob/AppData/Local/Temp/RtmpmW4XlM/response-families_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

## Bernoulli (accuracy)

Binary accuracy via a logit link: the intercept is the log-odds of a correct response, the
effect a log-odds difference between conditions.


``` r
acc <- demo("bernoulli", intercept = 0, effect = 0.5)
round(acc$by_group, 3)   # P(correct) by group
#>   control treatment 
#>     0.446     0.569
```

## Poisson (counts)

Counts via a log link (e.g. number of fixations, errors, or events). An intercept of 1.5
implies a base rate near `exp(1.5)`.


``` r
cts <- demo("poisson", intercept = 1.5, effect = 0.3)
round(cts$by_group, 2)   # mean count by group
#>   control treatment 
#>      3.85      5.18
table(cts$y)[1:8]
#> 
#>   0   1   2   3   4   5   6   7 
#>  41 225 466 698 711 659 466 334
```

## Ordinal (Likert)

Ordered categorical responses via a cumulative-logit model with user thresholds. The effect
shifts the latent distribution across the thresholds.


``` r
ord <- build_spec(list(
  name = "likert", seed = 1, design_kind = "between", n_subject = 4000,
  factor_name = "group", lev1 = "control", lev2 = "treatment",
  intercept = 0, effect = 0.8, family = "ordinal", resp_name = "rating",
  thresholds = "-2, -0.6, 0.6, 2"))
r <- simulate_design(ord)
round(prop.table(table(r$group, r$rating), 1), 2)   # category proportions by group
#>            
#>                1    2    3    4    5
#>   control   0.16 0.30 0.27 0.19 0.08
#>   treatment 0.09 0.18 0.29 0.28 0.17
```

## Beta (proportions)

Bounded proportions in (0, 1) via a mean–precision parameterization: the mean is
`logit⁻¹(η)` and `phi` is the precision (larger = tighter).


``` r
bt <- demo("beta", intercept = 0, effect = 0.8, phi = 8)
round(bt$by_group, 3)   # mean proportion by group
#>   control treatment 
#>     0.400     0.592
hist(bt$y, breaks = 40, col = "#2E8B57", border = "white", main = "Beta", xlab = "proportion")
```

![](C:/Users/pablob/AppData/Local/Temp/RtmpmW4XlM/response-families_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

## Choosing a family

| Outcome | Family | Scale of the effect |
|---|---|---|
| Continuous, symmetric | `gaussian` | identity |
| Reaction time | `shifted_lognormal` | log |
| Accuracy (0/1) | `bernoulli` | logit |
| Counts | `poisson` | log |
| Likert / ordered categories | `ordinal` | logit (cumulative) |
| Proportions in (0, 1) | `beta` | logit (mean) |

These match the families researchers fit in `lme4`, `glmmTMB`, and `brms`, so a design
simulated here lines up with the model you will fit. The point-and-click app exposes all six;
the full engine (continuous predictors, interactions, nesting, partial crossing) is available
by writing the specification directly.
