# Rescale a design to a target total variance

Adjust a specification so that the total variance of its linear
predictor, plus the residual variance of its response family, comes to
`target_var`. Calibrating to 1 puts the outcome on a unit scale, which
is what lets a region of practical equivalence or a smallest effect size
of interest be stated in standard-deviation units and read the same way
across designs. The total includes the variance the fixed effects
contribute, so a coefficient calibrated against it is in standard
deviations of the outcome and is not Cohen's d. With
`denominator = "random"`, the target applies to the total less that
fixed component, the variance that d divides by.

## Usage

``` r
calibrate_response(
  spec,
  target_var = 1,
  tune = c("sigma", "all"),
  scale_intercept = FALSE,
  denominator = c("total", "random")
)
```

## Arguments

- spec:

  A design specification (path or list).

- target_var:

  The variance to calibrate to, as `denominator` defines it. Defaults to
  1.

- tune:

  Either `"sigma"` (the default), which solves for the residual standard
  deviation alone, or `"all"`, which scales every variance-contributing
  term by a common factor.

- scale_intercept:

  With `tune = "all"`, whether the intercept is scaled with the rest.
  `FALSE` (the default) holds it, and `TRUE` reproduces pilotr 0.3.0 and
  0.3.1.

- denominator:

  Either `"total"` (the default), the total variance with the fixed
  component included, or `"random"`, the random-effect components plus
  the residual, as in Cohen's d.

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

`tune = "all"` holds the intercept and multiplies every coefficient and
every random-effect standard deviation, and the residual standard
deviation where there is one, by a common factor. For the families with
a free residual that leaves every ratio between components unchanged, so
it rescales the outcome without altering the design's character. It is
also the only option for the families whose residual is fixed by the
link.

For those, the residual cannot be rescaled at all, so the structural
part alone has to close the gap and the target has to exceed the link's
own contribution. A `bernoulli` or `ordinal` design carries a latent
residual of `pi^2 / 3`, about 3.29, so calibrating one to a total
variance of 1 is not merely difficult but impossible, and the error says
so rather than returning something plausible. For `poisson` and `beta`
the residual moves with the linear predictor, so the factor is solved
numerically; that costs one extra simulation rather than one per
candidate, because scaling the coefficients and the random-effect
standard deviations by `k` moves each row's linear predictor from `eta`
to exactly `b0 + k * (eta - b0)`, where `b0` is the intercept. The
search starts at `k = 0`, where every row sits at the intercept, so a
target at or below the residual variance that the intercept alone
implies is refused.

The intercept is held because it is a location. It adds nothing to the
fixed or random components, so scaling it changes none of them and only
moves the outcome. For `poisson` and `beta`, whose residual depends on
the location, scaling the intercept would reach the target partly by
moving the mean count or proportion. In the first example below, the
median reaction time of 600 ms stays near 600 ms. The intercept is the
location when the contrasts and predictors are centred, as codes of -0.5
and 0.5 are. Scaling the coefficient of a contrast coded 0 and 1, or of
a predictor whose mean is not zero, moves the location as well.

`scale_intercept = TRUE` multiplies the intercept by the same factor, as
pilotr 0.3.0 and 0.3.1 did. For `gaussian` and `exgaussian`, that is a
change of units, and no power or precision result for a coefficient
changes. Through a log or logit link, it changes the outcome itself. For
the lognormal families, it raises the response, less any shift, to the
power of the factor, which turned that median of 600 ms into 16.6 hours.
For `poisson`, it refuses targets that a held intercept reaches, and it
can return an intercept too large to simulate.

`denominator` says what `target_var` is a target for. `"total"`, the
default, is the total that
[`response_variance()`](https://pablobernabeu.github.io/pilotr/r/reference/response_variance.md)
reports, fixed component included. `"random"` is that total less the
fixed component, leaving the random-effect components and the residual.
In a balanced design, that is the variance within a condition. Westfall,
Kenny and Judd (2014) define Cohen's d for designs in which participants
respond to stimuli by dividing a difference between condition means by
the square root of that variance. Once it is calibrated to 1, the
coefficient of a two-level factor coded -0.5 and 0.5 is d, on the scale
of the linear predictor and conditional on any covariates. The fixed
component comes on top, so the total then exceeds `target_var`. With
`tune = "sigma"`, the residual standard deviation makes up what the
random components leave of the target.
[`response_variance()`](https://pablobernabeu.github.io/pilotr/r/reference/response_variance.md)
gives the conversion of a published d into a coefficient.

## References

Brysbaert, M. and Stevens, M. (2018). Power analysis and effect size in
mixed effects models: A tutorial. *Journal of Cognition*, 1(1), 9.
[doi:10.5334/joc.10](https://doi.org/10.5334/joc.10)

Westfall, J., Kenny, D. A. and Judd, C. M. (2014). Statistical power and
optimal design in experiments in which samples of participants respond
to samples of stimuli. *Journal of Experimental Psychology: General*,
143(5), 2020-2045.
[doi:10.1037/xge0000014](https://doi.org/10.1037/xge0000014)

## See also

[`response_variance()`](https://pablobernabeu.github.io/pilotr/r/reference/response_variance.md)
for the decomposition this works from.

## Examples

``` r
spec <- pilotr_example("crossed_mixed_rt")
calibrated <- calibrate_response(spec, target_var = 1, tune = "all")
round(response_variance(calibrated)$total, 6)
#> [1] 1
calibrated$fixed$intercept   # held at 6, so the median reaction time stays near 600 ms
#> [1] 6

# Calibrate the random components plus the residual, so that the coefficient of `cond`, coded
# -0.5 and 0.5, is Cohen's d on the log scale.
d_scale <- calibrate_response(spec, target_var = 1, tune = "all", denominator = "random")
d_scale$fixed$coefficients$cond
#> [1] 0.1498727
```
