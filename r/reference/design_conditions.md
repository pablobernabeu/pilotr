# Build a grid of fixed-effect coefficient sets

Produce the conditions needed to sweep an effect size with
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

  Optional named list of coefficients to include in every condition, the
  null condition among them.
  [`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md)
  already keeps the coefficients a condition does not name. `.base` is
  therefore needed only to hold a coefficient at a value other than the
  specification's, or for a condition used outside
  [`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md),
  such as one assigned to `fixed$coefficients` by hand.

## Value

A list of class `pilotr_conditions` holding one named list of
coefficients per condition. The class is how
[`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md)
recognises conditions to merge.

## Details

Named arguments give the values each effect should take, and are
recycled to a common length, so
`design_conditions(cond = c(0.02, 0.05), age = 0.1)` produces two
conditions, both with `age` at 0.1. The all-zero condition comes first
and is shared, since the Type I error rate is a property of the design
rather than of any one effect size.

[`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md)
merges each condition into the specification's own coefficients and
keeps their order, so any coefficient not named here keeps its value and
a sweep varies only the effects it names. It refuses a condition that
names a coefficient the specification lacks. Subsetting the result with
`[` keeps its class, and so does [`c()`](https://rdrr.io/r/base/c.html)
when the result comes first. The three conditions of
`c(design_conditions(cond = 0.02), design_conditions(age = 0.1, .null = FALSE))`
are therefore merged in the same way. A list assembled by other means
replaces the coefficients wholesale.

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
