# Rescale a design to a target total variance

Adjust a specification so that the total variance of its linear
predictor, plus the residual variance of its response family, comes to
`target_var`. Calibrating to 1 puts the outcome on a unit scale, which
is what lets a region of practical equivalence or a smallest effect size
of interest be stated in standard-deviation units and read the same way
across designs.

## Usage

``` r
calibrate_response(spec, target_var = 1, tune = c("sigma", "all"))
```

## Arguments

- spec:

  A design specification (path or list).

- target_var:

  The total variance to calibrate to. Defaults to 1.

- tune:

  Either `"sigma"` (the default), which solves for the residual standard
  deviation alone, or `"all"`, which scales every variance-contributing
  term by a common factor.

## Value

The specification, with the tuned parameters replaced.

## Details

`tune = "sigma"` holds the fixed effects and the random-effect standard
deviations where they are and solves for the residual standard
deviation. That keeps every effect size in the specification as written,
and is the right choice when those effects come from a pilot study or
from the literature. It fails when the structural variance already
exceeds the target, since no residual standard deviation can bring the
total down, and the error says so rather than returning a negative
variance.

`tune = "all"` multiplies the intercept, every coefficient and every
random-effect standard deviation, and the residual standard deviation
where there is one, by a common factor. For the families with a free
residual that leaves every ratio between components unchanged, so it
rescales the outcome without altering the design's character. It is also
the only option for the families whose residual is fixed by the link.

For those, the residual cannot be rescaled at all, so the structural
part alone has to close the gap and the target has to exceed the link's
own contribution. A `bernoulli` or `ordinal` design carries a latent
residual of `pi^2 / 3`, about 3.29, so calibrating one to a total
variance of 1 is not merely difficult but impossible, and the error says
so rather than returning something plausible. For `poisson` and `beta`
the residual moves with the linear predictor, so the factor is solved
numerically; that costs one extra simulation rather than one per
candidate, because scaling every term by `k` multiplies each row's
linear predictor by exactly `k`.

## See also

[`response_variance()`](https://pablobernabeu.github.io/pilotr/r/reference/response_variance.md)
for the decomposition this works from.

## Examples

``` r
spec <- pilotr_example("crossed_mixed_rt")
calibrated <- calibrate_response(spec, target_var = 1, tune = "all")
round(response_variance(calibrated)$total, 6)
#> [1] 1
```
