# Build the modelling data frame from a simulated data set and its specification

Add the analysis response column \`.y\` (log-transformed for the
lognormal families), the numeric contrast columns implied by the
categorical factors, and any interaction product columns (an \`a:b\`
coefficient becomes a column \`a_b\`).

## Usage

``` r
model_data(spec, d)
```

## Arguments

- spec:

  A design specification (path or list).

- d:

  A simulated data set, as returned by \[simulate_design()\].

## Value

The data frame \`d\` augmented with the \`.y\` response column and the
contrast and interaction columns required by the auto-derived model.

## Examples

``` r
spec <- build_spec(list(name = "d", seed = 1, design_kind = "between",
  factor_name = "group", lev1 = "a", lev2 = "b", n_subject = 20,
  intercept = 0, effect = 0.5, family = "gaussian", resp_name = "", sigma = 1))
head(model_data(spec, simulate_design(spec)))
#>   subject group   score      .y effect
#> 1       1     a  0.4240  0.4240   -0.5
#> 2       2     a -0.8216 -0.8216   -0.5
#> 3       3     a  0.8925  0.8925   -0.5
#> 4       4     a -1.1222 -1.1222   -0.5
#> 5       5     a -1.3471 -1.3471   -0.5
#> 6       6     a -1.0623 -1.0623   -0.5
```
