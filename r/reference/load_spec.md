# Load a design specification from a JSON file

Load a design specification from a JSON file

## Usage

``` r
load_spec(path, validate = TRUE)
```

## Arguments

- path:

  Path to a JSON design-specification file.

- validate:

  Whether to validate the specification after reading it. `TRUE` (the
  default) applies
  [`validate_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/validate_spec.md)
  with `strict = TRUE`; `FALSE` skips validation, and any other value is
  passed to
  [`validate_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/validate_spec.md)
  as its `strict` argument.

## Value

The specification as a nested list, with sub-lists left unsimplified so
that the structure round-trips exactly. Pass the result to
[`simulate_design()`](https://pablobernabeu.github.io/pilotr/r/reference/simulate_design.md).

## Details

The specification is validated by default, via
[`validate_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/validate_spec.md),
because several ways of getting one wrong produce plausible data and no
error at all: a mistyped coefficient key resolves to no column and so
silently sets that effect to zero, and a response parameter left over
from another family is ignored. Validation also refuses a specification
declaring a `spec_version` newer than this implementation understands,
and never reads such a file in part.

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
