# Simulate a data set from a design specification

Generate an analysis-ready data set from a portable design
specification: a linear predictor built from fixed effect sizes
(categorical contrasts, continuous predictors, and their interactions)
plus crossed by-subject and by-item random intercepts and slopes, mapped
through the chosen response family.

## Usage

``` r
simulate_design(spec, validate = TRUE)
```

## Arguments

- spec:

  A design specification, given either as a path to a JSON file or as an
  already-parsed list (for example from
  [`build_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/build_spec.md)
  or
  [`load_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/load_spec.md)).

- validate:

  Whether to validate the specification first, via
  [`validate_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/validate_spec.md).
  The default `TRUE` catches the errors that would otherwise pass
  silently, such as a mistyped coefficient key, which resolves to no
  column and so sets that effect to zero. Validation costs a few
  milliseconds, so the replicate loops behind the power and precision
  functions validate once and then pass `FALSE`; there is rarely a
  reason to set it directly.

## Value

A data frame with one row per observation, containing a `subject`
column, an optional `item` column, any grouping, factor, and
continuous-predictor columns, and the response column named by the
specification.

## Examples

``` r
spec <- build_spec(list(name = "demo", seed = 1, design_kind = "between",
  factor_name = "group", lev1 = "control", lev2 = "treatment", n_subject = 40,
  intercept = 100, effect = 5, family = "gaussian", resp_name = "", sigma = 10))
head(simulate_design(spec))
#>   subject   group    score
#> 1       1 control 104.2395
#> 2       2 control  91.7843
#> 3       3 control 108.9245
#> 4       4 control  88.7779
#> 5       5 control  86.5289
#> 6       6 control  89.3767
```
