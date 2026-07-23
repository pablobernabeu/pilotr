# Serialise a design specification to pretty-printed JSON

Serialise a design specification to pretty-printed JSON

## Usage

``` r
spec_json(spec)
```

## Arguments

- spec:

  A design specification (list), as produced by
  [`build_spec()`](https://pablobernabeu.github.io/pilotr/reference/build_spec.md).

## Value

A length-one character string containing the specification as
pretty-printed JSON, the portable artefact that the R and 'Python'
packages both consume.

## Examples

``` r
spec <- build_spec(list(name = "demo", seed = 1, design_kind = "between",
  factor_name = "group", lev1 = "a", lev2 = "b", n_subject = 20,
  intercept = 0, effect = 0.5, family = "gaussian", resp_name = "", sigma = 1))
cat(spec_json(spec))
#> {
#>   "name": "demo",
#>   "seed": 1,
#>   "units": {
#>     "subject": {
#>       "n": 20
#>     }
#>   },
#>   "factors": [
#>     {
#>       "name": "group",
#>       "levels": ["a", "b"],
#>       "contrasts": {
#>         "effect": [-0.5, 0.5]
#>       },
#>       "between": "subject"
#>     }
#>   ],
#>   "fixed": {
#>     "intercept": 0,
#>     "coefficients": {
#>       "effect": 0.5
#>     }
#>   },
#>   "response": {
#>     "family": "gaussian",
#>     "name": "score",
#>     "sigma": 1,
#>     "round": 4
#>   }
#> }
```
