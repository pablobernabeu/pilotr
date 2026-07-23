# Generate a self-contained, reproducible R script from a specification

Embed the specification as an R list literal (via `deparse`, which
round-trips exactly), so that the returned script reproduces the design
without any external file. This turns a design built in the no-code
application into a reproducible script; the application's Verify button
runs that script in a clean R session and confirms that it reproduces
the data bit-for-bit.

## Usage

``` r
generate_r_script(spec)
```

## Arguments

- spec:

  A design specification (list), as produced by
  [`build_spec()`](https://pablobernabeu.github.io/pilotr/reference/build_spec.md).

## Value

A length-one character string containing a runnable R script that loads
`pilotr`, embeds the specification, and simulates the data.

## Examples

``` r
spec <- build_spec(list(name = "demo", seed = 1, design_kind = "between",
  factor_name = "group", lev1 = "a", lev2 = "b", n_subject = 20,
  intercept = 0, effect = 0.5, family = "gaussian", resp_name = "", sigma = 1))
cat(generate_r_script(spec))
#> # Reproducible simulation exported by pilotr.
#> # install.packages("pilotr")   # once available; then run this script as-is.
#> library(pilotr)
#> 
#> spec <- list(name = "demo", seed = 1L, units = list(subject = list(n = 20L)), 
#>     factors = list(list(name = "group", levels = c("a", "b"), 
#>         contrasts = list(effect = c(-0.5, 0.5)), between = "subject")), 
#>     fixed = list(intercept = 0, coefficients = list(effect = 0.5)), 
#>     response = list(family = "gaussian", name = "score", sigma = 1, 
#>         round = 4L))
#> 
#> data <- simulate_design(spec)              # analysis-ready data frame
#> # write.csv(data, "data.csv", row.names = FALSE)
#> # pow  <- power_mixed(spec, n_sims = 200)   # simulation-based power + Type S/M
```
