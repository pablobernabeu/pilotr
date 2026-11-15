# Print a brms bridge

Writes the ready-to-fit `brms` model held in the object's `code`
element, which is the form the bridge is meant to be read in and the one
to copy into a script. The other elements (`formula`, `family`,
`priors`) are the same model in parts, for a caller assembling its own
code, and are left to [`str()`](https://rdrr.io/r/utils/str.html) or to
`$`.

## Usage

``` r
# S3 method for class 'pilotr_bridge'
print(x, ...)
```

## Arguments

- x:

  A `pilotr_bridge` object, as returned by
  [`brms_bridge()`](https://pablobernabeu.github.io/pilotr/r/reference/brms_bridge.md).

- ...:

  Ignored, present for consistency with the generic.

## Value

`x`, invisibly.

## Examples

``` r
spec <- build_spec(list(name = "d", seed = 1, design_kind = "between",
  factor_name = "g", lev1 = "a", lev2 = "b", n_subject = 20,
  intercept = 0, effect = 0.4, family = "gaussian", resp_name = "y", sigma = 1))
print(brms_bridge(spec))
#> library(brms)
#> fit <- brm(
#>   y ~ effect,
#>   data   = your_data,   # model_data(spec, simulate_design(spec)) adds the contrast columns
#>   family = gaussian(),
#>   prior  = c(
#>     prior(normal(0, 0.995), class = "b", coef = "effect")
#>   ),
#>   chains = 4, iter = 4000, warmup = 2000, cores = 4,
#>   control = list(adapt_delta = 0.95)
#> ) 
```
