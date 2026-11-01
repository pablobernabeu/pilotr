# Default response-column name for a family

Default response-column name for a family

## Usage

``` r
default_response_name(family)
```

## Arguments

- family:

  A response-family name, one of `"gaussian"`, `"lognormal"`,
  `"shifted_lognormal"`, `"exgaussian"`, `"bernoulli"`, `"poisson"`,
  `"ordinal"` or `"beta"`.

## Value

The conventional response-column name for that family (for example
`"RT"` for the reaction-time families `"lognormal"`,
`"shifted_lognormal"` and `"exgaussian"`), or `"outcome"` for an
unrecognised family.

## Examples

``` r
default_response_name("bernoulli")
#> [1] "accuracy"
```
