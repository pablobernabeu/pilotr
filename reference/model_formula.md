# Derive the lmer formula implied by a specification

Construct the mixed-model formula (with response `.y`, as produced by
[`model_data()`](https://pablobernabeu.github.io/pilotr/reference/model_data.md))
from the fixed-effect coefficients and random-effects structure of a
specification. Interaction coefficients written `a:b` become formula
terms `a_b`.

## Usage

``` r
model_formula(spec)
```

## Arguments

- spec:

  A design specification (path or list).

## Value

A [stats::formula](https://rdrr.io/r/stats/formula.html) object suitable
for fitting with `lme4` or `lmerTest`.

## Examples

``` r
spec <- build_spec(list(name = "d", seed = 1, design_kind = "within",
  include_items = TRUE, n_subject = 10, n_item = 8, factor_name = "cond",
  lev1 = "a", lev2 = "b", intercept = 6, effect = 0.05,
  subj_int_sd = 0.12, subj_slope_sd = 0.04, subj_corr = 0.2,
  item_int_sd = 0.08, item_slope_sd = 0.02, item_corr = -0.1,
  family = "shifted_lognormal", resp_name = "", sigma = 0.3, shift = 200))
model_formula(spec)
#> .y ~ effect + (1 + effect | subject) + (1 + effect | item)
```
