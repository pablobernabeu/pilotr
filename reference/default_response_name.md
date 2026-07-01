# Default response-column name for a family

Default response-column name for a family

## Usage

``` r
default_response_name(family)
```

## Arguments

- family:

  A response-family name, one of \`"gaussian"\`,
  \`"shifted_lognormal"\`, \`"bernoulli"\`, \`"poisson"\`,
  \`"ordinal"\`, or \`"beta"\`.

## Value

The conventional response-column name for that family (for example
\`"RT"\` for \`"shifted_lognormal"\`), or \`"outcome"\` for an
unrecognised family.

## Examples

``` r
default_response_name("bernoulli")
#> [1] "accuracy"
```
