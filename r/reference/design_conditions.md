# Build a grid of fixed-effect coefficient sets

Produce the list of `fixed$coefficients` objects needed to sweep an
effect size with
[`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md),
including a condition in which every named effect is zero, so that the
same run shows both what a design detects and how often it declares
something when there is nothing to find.

## Usage

``` r
design_conditions(..., .null = TRUE, .base = NULL)
```

## Arguments

- ...:

  Named numeric vectors, one per coefficient to vary.

- .null:

  Whether to prepend a condition with every named effect set to zero.
  `TRUE` by default.

- .base:

  Optional named list of coefficients to hold fixed in every condition,
  for effects the sweep does not vary.

## Value

A list of named lists, each suitable as a `fixed$coefficients` value.

## Details

Named arguments give the values each effect should take, and are
recycled to a common length, so
`design_conditions(cond = c(0.02, 0.05), age = 0.1)` produces two
conditions, both with `age` at 0.1. The all-zero condition comes first
and is shared, since the Type I error rate is a property of the design
rather than of any one effect size.

Any coefficient the specification has that is not named here is left at
its own value, so a sweep varies only the effects it names.

## See also

[`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md),
which consumes this.

## Examples

``` r
design_conditions(effect = c(0.03, 0.06))
#> [[1]]
#> [[1]]$effect
#> [1] 0
#> 
#> 
#> [[2]]
#> [[2]]$effect
#> [1] 0.03
#> 
#> 
#> [[3]]
#> [[3]]$effect
#> [1] 0.06
#> 
#> 
design_conditions(cond = c(0.02, 0.05), age = 0.1, .null = FALSE)
#> [[1]]
#> [[1]]$cond
#> [1] 0.02
#> 
#> [[1]]$age
#> [1] 0.1
#> 
#> 
#> [[2]]
#> [[2]]$cond
#> [1] 0.05
#> 
#> [[2]]$age
#> [1] 0.1
#> 
#> 
```
