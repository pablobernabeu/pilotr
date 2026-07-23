# Load a design specification from a JSON file

Load a design specification from a JSON file

## Usage

``` r
load_spec(path)
```

## Arguments

- path:

  Path to a JSON design-specification file.

## Value

The specification as a nested list, with sub-lists left unsimplified so
that the structure round-trips exactly. Pass the result to
[`simulate_design()`](https://pablobernabeu.github.io/pilotr/reference/simulate_design.md).

## Examples

``` r
spec <- build_spec(list(name = "demo", seed = 1, design_kind = "between",
  factor_name = "group", lev1 = "a", lev2 = "b", n_subject = 20,
  intercept = 0, effect = 0.5, family = "gaussian", resp_name = "", sigma = 1))
f <- tempfile(fileext = ".json")
writeLines(spec_json(spec), f)
identical(simulate_design(load_spec(f)), simulate_design(spec))
#> [1] TRUE
```
