# Print a simulation-based power result

Shows the fit accounting, with the share of the replicates analysed at
each rung when the result comes from a ladder, and then two tables with
one row per focal effect. The first gives each power estimate beside its
Monte Carlo standard error and Wilson interval, so that the precision of
the estimate is as visible as the estimate. The second gives Type S
beside its Wilson interval and Type M beside its Monte Carlo standard
error, both taken over the significant replicates alone.

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

  Number of decimal places for the reported rates and for the interval
  bounds.

- ...:

  Ignored, present for consistency with the generic.

## Value

`x`, invisibly.

## Details

A result whose `focal` named sets of coefficients lists the coefficients
of each set before the first table, which has one row per set and no
true value. It has no second table, since a joint test has no single
estimate from which Type S or Type M could be taken.
