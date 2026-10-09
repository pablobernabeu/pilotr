# Print a simulation-based power result

Shows the fit accounting and then one row per focal effect, with each
power estimate beside its Monte Carlo standard error and Wilson
interval, so that the precision of the estimate is as visible as the
estimate.

## Usage

``` r
# S3 method for class 'pilotr_power'
print(x, digits = 3, ...)
```

## Arguments

- x:

  A `pilotr_power` object, as returned by
  [`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md).

- digits:

  Number of significant digits for the reported rates.

- ...:

  Ignored, present for consistency with the generic.

## Value

`x`, invisibly.
