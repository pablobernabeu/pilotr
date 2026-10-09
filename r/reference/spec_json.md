# Serialise a design specification to pretty-printed JSON

Serialise a design specification to pretty-printed JSON

## Usage

``` r
spec_json(spec)
```

## Arguments

- spec:

  A design specification (list), as produced by
  [`build_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/build_spec.md).

## Value

A length-one character string containing the specification as
pretty-printed JSON, the portable artefact that the R and 'Python'
packages both consume.

## Details

Each number is written at the shortest width that reads back as exactly
the same IEEE-754 double, which is 15 significant digits for most values
a user types and at most 17 for the rest. The JSON file is the portable
artefact that the R and 'Python' implementations both read, so a fixed
lower precision makes the specification itself a source of
cross-language divergence: at 15 digits throughout, a coefficient of
`1/3` came back as `0.33333333333333298`, and over a sample of 214
doubles 189 failed to round-trip. Writing every number at 17 instead
would round-trip but read badly, turning an effect typed as 0.3 into
`0.29999999999999999`.

Preferring the shorter width is not only a matter of appearance. A build
of R without long-double arithmetic reads 17 significant digits
inexactly, because the digits overflow the 53-bit mantissa before the
decimal exponent is applied, so the few values that need that width are
the ones such a build cannot recover. Anything shorter it reads
correctly.

## Examples

``` r
spec <- build_spec(list(name = "demo", seed = 1, design_kind = "between",
  factor_name = "group", lev1 = "a", lev2 = "b", n_subject = 20,
  intercept = 0, effect = 0.5, family = "gaussian", resp_name = "", sigma = 1))
cat(spec_json(spec))
#> {
#>   "spec_version": "0.3",
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
