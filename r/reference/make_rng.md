# Create a shared cross-language random-number generator

Build the combined linear congruential generator (L'Ecuyer 1988) used by
both the R and Python implementations, so that a given specification and
seed produce identical data in either language. The draw-order contract
is documented in the specification at
<https://github.com/pablobernabeu/pilotr/blob/main/spec/SPEC.md>.

## Usage

``` r
make_rng(seed)
```

## Arguments

- seed:

  A single number used to seed the generator (coerced to a non-negative
  integer).

## Value

A list of three functions: `uniform()` returns one standard-uniform
draw, `normal()` returns one standard-normal draw, and `normals(k)`
returns a length-`k` numeric vector of standard-normal draws.

## Examples

``` r
rng <- make_rng(1)
rng$uniform()
#> [1] 0.7498299
rng$normals(3)
#> [1] -0.5715727  1.1424531 -0.8722080
```
