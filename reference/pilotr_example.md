# Example design specifications shipped with pilotr

pilotr ships one ready-to-run specification per design family, as JSON,
in the package's `examples/` directory. These are the same files that
drive the Python twin and the no-code app, so a design authored in one
place runs unchanged in the others. `pilotr_example()` lists them, or
returns the path to one for
[`load_spec()`](https://pablobernabeu.github.io/pilotr/reference/load_spec.md).

## Usage

``` r
pilotr_example(name = NULL)
```

## Arguments

- name:

  The base name of an example, with or without the `.json` extension,
  for example `"between_2group_gaussian"`. When `NULL` (the default),
  the available example names are returned instead of a path.

## Value

When `name` is `NULL`, a character vector of the available example
names. Otherwise, the full path to that example's JSON file, ready to
pass to
[`load_spec()`](https://pablobernabeu.github.io/pilotr/reference/load_spec.md).

## See also

[`load_spec()`](https://pablobernabeu.github.io/pilotr/reference/load_spec.md)
to read a specification and
[`simulate_design()`](https://pablobernabeu.github.io/pilotr/reference/simulate_design.md)
to simulate from it.

## Examples

``` r
pilotr_example()                       # the available examples
#> [1] "beta_proportion"         "between_2group_gaussian"
#> [3] "crossed_mixed_rt"        "nested_clusters"        
#> [5] "ordinal_likert_between"  "partial_crossing"       
#> [7] "poisson_counts_between"  "reading_time_continuous"
spec <- load_spec(pilotr_example("between_2group_gaussian"))
head(simulate_design(spec))
#>   subject   group    score
#> 1       1 control  95.7157
#> 2       2 control  90.0921
#> 3       3 control 119.1958
#> 4       4 control  86.4388
#> 5       5 control 103.4220
#> 6       6 control  92.0864
```
