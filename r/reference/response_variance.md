# Variance components of the linear predictor

Decompose the variance of a design's linear predictor into the part
contributed by the fixed effects, the part contributed by each grouping
factor's random effects, and the residual variance added by the response
family. Useful for putting a region of practical equivalence or a
smallest effect size of interest on a known scale, and for seeing which
term dominates a design before committing to it.

## Usage

``` r
response_variance(spec)
```

## Arguments

- spec:

  A design specification (path or list).

## Value

A named list of variance components: `fixed`, one entry per grouping
factor, `residual`, and `total` (their sum).

## Details

The fixed component is the variance, across rows, of the linear
predictor with every random-effect standard deviation set to zero. It is
read off the design the specification actually produces, so it needs no
assumption about how a between-unit factor divides the units or about
whether the predictors are independent.

Each grouping factor's component is the average over rows of
`x' Sigma x`, where `Sigma` is the covariance the specification asks for
and `x` collects that row's values of the intercept and each
random-slope column. This is exact for the realised design: it averages
over the random-effect distribution analytically, which matters because
drawing from it and estimating a variance from the drawn effects of,
say, 30 subjects carries a sampling error of around a quarter of the
component itself, far too much to calibrate a region of practical
equivalence against.

The cost is one simulation for the fixed part plus one per distinct
random-slope column, which for a typical design is a handful.

Every component is on the scale the linear predictor lives on, which is
what makes them comparable. For `lognormal` and `shifted_lognormal` that
is the log of the response, which is also the scale the auto-derived
analysis model works on. For `bernoulli`, `poisson`, `ordinal` and
`beta` it is the latent scale behind the link, so the residual is the
distribution-specific variance used to compute the intraclass
correlation and R-squared of a generalised mixed model (Nakagawa,
Johnson and Schielzeth, 2017).

Three of those four are derived from the process pilotr simulates and
are exact for it. A `bernoulli` row is drawn as `1[u < invlogit(eta)]`,
equivalently `logit(u) < eta`, so the latent error is a standard
logistic variate of variance `pi^2 / 3`, and `ordinal` compares the same
uniform against cumulative thresholds and inherits it. For `beta`,
`logit(Y)` has variance `trigamma(a) + trigamma(b)` exactly.

`poisson` is the one approximation, because a count of zero has no
logarithm and a log-scale variance cannot be measured directly. The
trigamma form is used, following the recommendation for the log link. It
agrees closely with the alternatives once the mean exceeds about 5 and
diverges sharply below it, where a log-scale variance is barely
meaningful: at a mean of 0.5 it gives 4.93 against 1.10 for the
lognormal approximation. Read the Poisson residual as an order of
magnitude when counts are rare.

The total includes the variance the fixed effects contribute. It is the
denominator of an R-squared, and of an effect or a region of practical
equivalence stated in standard deviations of the outcome, but it is not
the denominator of Cohen's d. For designs in which participants respond
to stimuli, Westfall, Kenny and Judd (2014) define d as the difference
between two condition means divided by the standard deviation within a
condition. The square of that standard deviation sums the random
intercept variances, the residual and the random slope variances, each
slope variance weighted by the squared contrast code. The fixed effects
are left out (see also Brysbaert and Stevens, 2018). In a balanced
design with codes that sum to zero, that is `total - fixed`. Covariates
belong to the fixed part, so a d defined this way is conditional on
them.

A published d therefore converts to a coefficient as
`b = d * sqrt(total - fixed) / (c2 - c1)`, where `c2 - c1` is the
difference between the two contrast codes. For codes of -0.5 and 0.5,
that difference is 1. The d is on the scale of the linear predictor, so
for the families with a link it is a latent-scale d.
[`calibrate_response()`](https://pablobernabeu.github.io/pilotr/r/reference/calibrate_response.md)
with `denominator = "random"` rescales a design so that `total - fixed`
is 1, after which the coefficient of a contrast coded -0.5 and 0.5 is d.

## References

Brysbaert, M. and Stevens, M. (2018). Power analysis and effect size in
mixed effects models: A tutorial. *Journal of Cognition*, 1(1), 9.
[doi:10.5334/joc.10](https://doi.org/10.5334/joc.10)

Nakagawa, S., Johnson, P. C. D. and Schielzeth, H. (2017). The
coefficient of determination R2 and intra-class correlation coefficient
from generalized linear mixed-effects models revisited and expanded.
*Journal of the Royal Society Interface*, 14(134), 20170213.
[doi:10.1098/rsif.2017.0213](https://doi.org/10.1098/rsif.2017.0213)

Westfall, J., Kenny, D. A. and Judd, C. M. (2014). Statistical power and
optimal design in experiments in which samples of participants respond
to samples of stimuli. *Journal of Experimental Psychology: General*,
143(5), 2020-2045.
[doi:10.1037/xge0000014](https://doi.org/10.1037/xge0000014)

## See also

[`calibrate_response()`](https://pablobernabeu.github.io/pilotr/r/reference/calibrate_response.md)
to rescale a design to a target total variance.

## Examples

``` r
spec <- pilotr_example("crossed_mixed_rt")
response_variance(spec)
#> $fixed
#> [1] 0.0006254343
#> 
#> $subject
#> [1] 0.0148
#> 
#> $item
#> [1] 0.0065
#> 
#> $residual
#> [1] 0.09
#> 
#> $total
#> [1] 0.1119254
#> 

# The coefficient that gives the related-unrelated contrast, coded -0.5 and 0.5, a published d
# of 0.2 on the log scale of this design.
v <- response_variance(spec)
0.2 * sqrt(v$total - v$fixed) / (0.5 - (-0.5))
#> [1] 0.06672331

# Every family reports a residual, including those whose outcome is discrete, so the components
# are a complete decomposition and their ratios read as the design's intraclass correlations.
v <- response_variance(pilotr_example("ordinal_likert_between"))
round(1 - v$residual / v$total, 3)   # share of latent variance that is structural
#> [1] 0.071
```
