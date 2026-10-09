# Precision and ROPE design analysis at a fixed sample size

A fast frequentist analogue of a Bayesian
highest-density-interval-versus-ROPE design analysis. Across Monte Carlo
replicates, fit the model and record, for each focal fixed effect,
whether its 95% confidence interval falls entirely outside a region of
practical equivalence (a practically meaningful effect) or entirely
inside it (practical equivalence to zero), along with the expected
interval width. Requires the `lme4` package.

## Usage

``` r
precision_design(
  spec,
  focal = NULL,
  formula = NULL,
  prep = NULL,
  rope = 0.05,
  n_sims = 100,
  workers = 1
)
```

## Arguments

- spec:

  A design specification (path or list).

- focal:

  The focal effects. `NULL`, the default, analyses every coefficient in
  the specification and takes the true values from it. A named numeric
  vector maps coefficient names to their true values, and a character
  vector names them without their true values. Interaction effects
  follow the model's column naming, so a specification key `a:b` is the
  focal name `a_b`.

- formula:

  Optional `lme4` formula; if `NULL` it is derived from the
  specification via
  [`model_formula()`](https://pablobernabeu.github.io/pilotr/r/reference/model_formula.md),
  and the response family chooses the fitter. A formula given here is
  fitted as a linear model whatever the family (see Details).

- prep:

  Optional function mapping a simulated data set to the modelling data;
  if `NULL` it is derived via
  [`model_data()`](https://pablobernabeu.github.io/pilotr/r/reference/model_data.md),
  which log-transforms the outcome and builds the contrast and
  interaction columns, so focal names follow the auto-formula
  (interactions written as `a_b`).

- rope:

  Half-width of the region of practical equivalence; an effect with
  `abs(beta) < rope` is treated as practically equivalent to zero. Set
  it clearly narrower than the smallest effect worth detecting, because
  the probability of a determinate meaningful decision about an effect
  no larger than `rope` cannot rise above 0.5 however large the sample.

- n_sims:

  Number of Monte Carlo replicates. `p_meaningful` and `p_equivalent`
  are proportions over the replicates that returned an estimate
  (`n_returned`), so they carry a Monte Carlo standard error of about
  `sqrt(p * (1 - p) / n_sims)` and move in coarse steps when `n_sims` is
  small. At least 200 replicates are advisable for real planning.

- workers:

  Number of local worker processes over which to spread the replicates.
  The default of 1 runs serially. Because every replicate takes its own
  seed from
  [`replicate_seeds()`](https://pablobernabeu.github.io/pilotr/r/reference/replicate_seeds.md),
  any worker count returns results identical to a serial run.

## Value

A data frame with one row per focal effect and columns `param`, `true`,
`mean_ci_width`, `p_meaningful`, `p_meaningful_mcse`, `p_meaningful_lo`,
`p_meaningful_hi`, `p_equivalent`, `p_equivalent_mcse`,
`p_equivalent_lo`, `p_equivalent_hi`, `n_attempted`, `n_returned`,
`n_converged`, `n_singular`, `n_warning` and `fitter`. Each decision
proportion is reported with its Monte Carlo standard error (`*_mcse`)
and Wilson interval bounds (`*_lo`, `*_hi`), because a proportion over a
finite number of replicates is an estimate rather than a fact. The
interval behind `mean_ci_width` and the ROPE decisions is the Wald
approximation described in Details. `fitter` names the fitter and the
interval, such as `"lme4::glmer (binomial), Wald z"`. For an `ordinal`
or `beta` response fitted on its own scale, the decision columns are
`NA` (see Details).

The decision proportions are taken over `n_returned`, the replicates
that produced an estimate. The remaining counts separate the fit
outcomes, because a fit can return a usable estimate while still being
boundary-singular or carrying a convergence warning: `n_converged`
counts replicates with neither, `n_singular` those where
[`lme4::isSingular()`](https://rdrr.io/pkg/lme4/man/isSingular.html) was
true, and `n_warning` those with a warning or optimiser convergence
message other than lme4's singular-fit notice, which `n_singular`
already counts. Singular and warning fits are retained, since their
fixed-effect estimates remain interpretable and discarding them would
bias the result. A large `n_singular` means the model being fitted is
richer than the design can support at that sample size, which is common
in crossed designs (Bates et al., 2015; Matuschek et al., 2017).

## Details

The interval is a Wald approximation: the estimate plus or minus 1.96
standard errors from the model's variance-covariance matrix. This
fixed-z interval is chosen for speed and for comparability across
replicates; in small samples it is somewhat narrower than a
Satterthwaite t interval, so `p_meaningful` and `mean_ci_width` are
slightly optimistic at small sample sizes.

The fitter follows the response family, as in
[`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md).
Since no p-value is needed, a `gaussian`, `lognormal`,
`shifted_lognormal` or `exgaussian` response is fitted by
[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html). A `bernoulli`
or `poisson` response is fitted by
[`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) with the
binomial or Poisson family, and a model with no random terms by
[`stats::lm()`](https://rdrr.io/r/stats/lm.html) or
[`stats::glm()`](https://rdrr.io/r/stats/glm.html). The interval is
therefore on the scale of the specification's coefficients, and so is
the region of practical equivalence that `rope` sets. That is the log
scale for the lognormal families and `poisson`, and the logit scale for
`bernoulli`. A `glmer()` fit has the Wald interval that Bolker et al.
(2009) recommend for a GLMM without overdispersion. With few subjects or
items, it is likewise narrower than it should be, since the test it
inverts rejects too often when the clusters are few (Li and Redden,
2015). The `fitter` column names the fitter.

pilotr has no frequentist model for an `ordinal` or `beta` response yet.
Such a response is fitted by the linear model on its own scale, where a
region of practical equivalence set on the logit scale of the
coefficients means nothing. The function warns, and withholds the
decision probabilities with their Monte Carlo standard errors and bounds
(`NA`). `mean_ci_width` stays, as the width of an interval on the
response scale.
[`generate_design_analysis()`](https://pablobernabeu.github.io/pilotr/r/reference/generate_design_analysis.md)
applies an interval-and-ROPE rule on the link scale, in a Bayesian
model.

A `formula` given here is fitted as written, by
[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html), or by
[`stats::lm()`](https://rdrr.io/r/stats/lm.html) when it has no random
terms, whatever the family. With such a formula, or a `prep` that
changes the response, and `focal = NULL`, the true values still come
from the specification and are on its scale, which need not be the
model's. A named numeric `focal` gives the true values on the model's
scale.

## References

Bates, D., Kliegl, R., Vasishth, S. and Baayen, H. (2015). Parsimonious
mixed models. *arXiv*.
[doi:10.48550/arXiv.1506.04967](https://doi.org/10.48550/arXiv.1506.04967)

Matuschek, H., Kliegl, R., Vasishth, S., Baayen, H. and Bates, D.
(2017). Balancing Type I error and power in linear mixed models.
*Journal of Memory and Language*, 94, 305-315.
[doi:10.1016/j.jml.2017.01.001](https://doi.org/10.1016/j.jml.2017.01.001)

Bolker, B. M., Brooks, M. E., Clark, C. J., Geange, S. W., Poulsen, J.
R., Stevens, M. H. H. and White, J.-S. S. (2009). Generalized linear
mixed models: A practical guide for ecology and evolution. *Trends in
Ecology & Evolution*, 24(3), 127-135.
[doi:10.1016/j.tree.2008.10.008](https://doi.org/10.1016/j.tree.2008.10.008)

Li, P. and Redden, D. T. (2015). Comparing denominator degrees of
freedom approximations for the generalized linear mixed model in
analyzing binary outcome in small sample cluster-randomized trials. *BMC
Medical Research Methodology*, 15, 38.
[doi:10.1186/s12874-015-0026-x](https://doi.org/10.1186/s12874-015-0026-x)

## Examples

``` r
# \donttest{
if (requireNamespace("lme4", quietly = TRUE)) {
  spec <- build_spec(list(name = "pr", seed = 1, design_kind = "within",
    include_items = TRUE, n_subject = 12, n_item = 12, factor_name = "cond",
    lev1 = "a", lev2 = "b", intercept = 6, effect = 0.05,
    subj_int_sd = 0.12, subj_slope_sd = 0.04, subj_corr = 0.2,
    item_int_sd = 0.08, item_slope_sd = 0.02, item_corr = -0.1,
    family = "shifted_lognormal", resp_name = "", sigma = 0.3, shift = 200))
  # n_sims is small so the example runs quickly. Use 200 or more for real planning.
  precision_design(spec, focal = c(effect = 0.05), rope = 0.02, n_sims = 10)
}
#>    param true mean_ci_width p_meaningful p_meaningful_mcse p_meaningful_lo
#> 1 effect 0.05     0.1526052          0.1        0.09486833      0.01787621
#>   p_meaningful_hi p_equivalent p_equivalent_mcse p_equivalent_lo
#> 1         0.40415            0                 0               0
#>   p_equivalent_hi n_attempted n_returned n_converged n_singular n_warning
#> 1       0.2775328          10         10           0         10         0
#>               fitter
#> 1 lme4::lmer, Wald z
# }
```
