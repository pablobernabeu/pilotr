# Derive a brms formula, family, and priors from a design spec

Derive a brms formula, family, and priors from a design spec

## Usage

``` r
brms_bridge(spec, prior_scale = 0.5, interaction_scale = NULL)
```

## Arguments

- spec:

  a design spec (path or list).

- prior_scale:

  SD of the Normal prior on fixed main effects (standardised scale).

- interaction_scale:

  SD of the Normal prior on interaction terms (default prior_scale/2).

## Value

Invisibly, a list with elements `formula`, `family`, `priors`, and
`code`; the `code` element (a ready-to-fit `brms` model) is also printed
to the console.

## Examples

``` r
spec <- build_spec(list(name = "d", seed = 1, design_kind = "within",
  include_items = TRUE, n_subject = 20, n_item = 12, factor_name = "cond",
  lev1 = "a", lev2 = "b", intercept = 6, effect = 0.05,
  subj_int_sd = 0.12, subj_slope_sd = 0.04, subj_corr = 0.2,
  item_int_sd = 0.08, item_slope_sd = 0.02, item_corr = -0.1,
  family = "shifted_lognormal", resp_name = "RT", sigma = 0.3, shift = 200))
bridge <- brms_bridge(spec)
#> library(brms)
#> fit <- brm(
#>   RT ~ effect + (1 + effect | subject) + (1 + effect | item),
#>   data   = your_data,
#>   family = shifted_lognormal(),
#>   prior  = c(
#>     prior(normal(0, 2.5), class = "Intercept"),
#>     prior(normal(0, 0.5), class = "b", coef = "effect"),
#>     prior(normal(0, 1), class = "sd"),
#>     prior(lkj(2), class = "cor")
#>   ),
#>   chains = 4, iter = 4000, warmup = 2000, cores = 4,
#>   control = list(adapt_delta = 0.95)
#> ) 
bridge$formula
#> [1] "RT ~ effect + (1 + effect | subject) + (1 + effect | item)"
```
