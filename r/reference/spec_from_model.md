# Build a design specification from a fitted mixed model

Read a design specification off a linear mixed model already fitted with
[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) or
[`lmerTest::lmer()`](https://rdrr.io/pkg/lmerTest/man/lmer.html). The
fixed effects, the random-effect standard deviations and correlations,
and the residual standard deviation are taken from the fit, and the
numbers of subjects and items may be raised at the same time, so that a
pilot study or a published model becomes the starting point of a power
analysis rather than a set of numbers to invent. Requires the `lme4`
package.

## Usage

``` r
spec_from_model(
  fit,
  name = NULL,
  seed = 1,
  n_subject = NULL,
  n_item = NULL,
  family = NULL,
  round = NULL
)
```

## Arguments

- fit:

  A fitted linear mixed model, of class `lmerMod` (from
  [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html)) or
  `lmerModLmerTest` (from
  [`lmerTest::lmer()`](https://rdrr.io/pkg/lmerTest/man/lmer.html)).
  Anything else is refused with a message saying what to pass instead.

- name:

  A label for the returned specification. Defaults to `"from_model"`.

- seed:

  The master seed of the returned specification. Defaults to 1.

- n_subject:

  Number of subjects for the returned specification. Defaults to the
  number of levels the fit actually had, and is normally raised above
  it, since scaling a pilot design up is the point of reading a
  specification off a pilot fit.

- n_item:

  Number of items, treated the same way as `n_subject`. An error when
  the model has no second crossed grouping factor to act as items,
  because there is then no item unit to resize.

- family:

  The response family of the returned specification. `NULL`, the
  default, gives `"gaussian"`, which is what a `lmer` fit is on the
  scale it was fitted on. A single string names another family, and a
  list such as `list(family = "shifted_lognormal", shift = 200)` also
  supplies the parameters of that family which a Gaussian fit cannot
  provide.

- round:

  Decimal places for the simulated response, passed through to
  `response.round`. `NULL`, the default, leaves the response unrounded.

## Value

A design specification as a nested list, carrying `spec_version` and
validated with
[`validate_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/validate_spec.md),
so it can be passed directly to
[`simulate_design()`](https://pablobernabeu.github.io/pilotr/r/reference/simulate_design.md),
[`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md)
or
[`spec_json()`](https://pablobernabeu.github.io/pilotr/r/reference/spec_json.md).
Two attributes record the readings that the fit did not settle on its
own: `group_mapping`, a named character vector giving the specification
unit each of the fit's grouping factors became, and `column_kinds`, a
named character vector giving each model-frame column the classification
`"factor"` or `"predictor"`. Both are attributes rather than fields so
that the specification itself stays within the portable schema.

## Details

Settling on plausible random-effect standard deviations with nothing to
read them from is the step that most often stops a simulation-based
power analysis before it starts, and the tutorials on the method
converge on the same remedy, which is to take the variance components
from a pilot fit or from a published model (Green and MacLeod, 2016;
Kumle et al., 2021; DeBruine and Barr, 2021). This function performs
that transfer, and the specification it returns can be enlarged and
passed straight to
[`simulate_design()`](https://pablobernabeu.github.io/pilotr/r/reference/simulate_design.md)
or
[`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md).

Every number is taken on the scale the model was fitted on, and the
returned family is Gaussian by default, because that is what a `lmer`
fit is. A model of log reaction times therefore yields an intercept,
coefficients and residual standard deviation on the log scale, which is
the right thing for the `lognormal` and `shifted_lognormal` families,
whose linear predictor lives on that scale as well. Ask for one of those
with `family`, in the list form when the family needs a parameter the
fit cannot supply, as in
`family = list(family = "shifted_lognormal", shift = 200)`.

Whether a numeric column of the model frame was a contrast-coded factor
or a continuous covariate is not recorded anywhere in a fitted model, so
it is inferred here. A column that is a factor, a character vector or a
logical vector is treated as categorical and its own contrast coding is
read with
[`stats::contrasts()`](https://rdrr.io/r/stats/contrasts.html), which
makes the emitted contrast-column names agree with the coefficient names
lme4 produced. A numeric column with exactly two distinct values is also
treated as categorical, on the reasoning that a two-valued numeric
predictor in a factorial experiment is a coded factor far more often
than it is a covariate. Any other numeric column becomes a continuous
predictor, with `mean` and `sd` taken from the data, one value per unit
for a unit-level variable so that an unbalanced design does not distort
them. This is a heuristic and nothing more, it is reported by a message
on every call, and a genuine two-valued covariate has to be moved from
`factors` to `predictors` by hand.

A column is placed as `between` a unit when it holds one value within
every member of that unit, and as `vary_within` when it takes several
values inside every unit; the same test gives a continuous predictor its
`varies_by`, which becomes `"observation"` when the predictor varies
inside every unit. A model that codes one factor twice, once as a factor
for its fixed effect and once as a numeric contrast for a random slope,
gives two separate specification terms, because the two columns are
separate columns in the model frame and nothing in the fit ties them
together.

Interactions need one further step. lme4 writes the interaction of two
model-frame columns as `a:b`, which is already the specification's
convention, but
[`model_data()`](https://pablobernabeu.github.io/pilotr/r/reference/model_data.md)
gives an interaction its own product column named `a_b`, so a
specification built from a fit of pilotr's own modelling data would
otherwise acquire a spurious independent term. Such a column is
recognised by checking the product identity in the data rather than by
reading its name, which both avoids mistaking a column called `z_freq`
for an interaction and settles where to split a name with several
underscores. A recognised product column is reported and re-keyed to
`a:b`, and it contributes no factor or predictor of its own.

Grouping factors named `subject` and `item` are used as they stand.
Otherwise the one with the most levels becomes `subject`, the largest
remaining factor that genuinely crosses it becomes `item`, and every
other grouping factor becomes an extra `random` entry with the `over`
and `n` fields that pilotr's additional grouping factors take, `over`
being decided by which unit the factor partitions. Any renaming is
reported by a message and recorded in the `group_mapping` attribute of
the result. When each subject saw only some of the items, the item unit
gains a `per_subject` count, so that the recovered design keeps the
partial crossing of the original rather than silently becoming fully
crossed.

A boundary-singular pilot fit is carried across as it stands, which
means a variance estimated at zero or a correlation estimated at exactly
plus or minus one. Those are faithful readings of the fit rather than
defects, but they are also the sign that the random-effect structure was
richer than the pilot could support, and a power analysis resting on
them will inherit that (Bates et al., 2015). Widening the design, or
simplifying the structure before refitting, is the remedy.

## References

Bates, D., Kliegl, R., Vasishth, S. and Baayen, H. (2015). Parsimonious
mixed models. *arXiv*.
[doi:10.48550/arXiv.1506.04967](https://doi.org/10.48550/arXiv.1506.04967)

DeBruine, L. M. and Barr, D. J. (2021). Understanding mixed-effects
models through data simulation. *Advances in Methods and Practices in
Psychological Science*, 4(1).
[doi:10.1177/2515245920965119](https://doi.org/10.1177/2515245920965119)

Green, P. and MacLeod, C. J. (2016). SIMR: an R package for power
analysis of generalized linear mixed models by simulation. *Methods in
Ecology and Evolution*, 7(4), 493-498.
[doi:10.1111/2041-210X.12504](https://doi.org/10.1111/2041-210X.12504)

Kumle, L., Vo, M. L.-H. and Draschkow, D. (2021). Estimating power in
(generalized) linear mixed models: An open introduction and tutorial in
R. *Behavior Research Methods*, 53, 2528-2543.
[doi:10.3758/s13428-021-01546-0](https://doi.org/10.3758/s13428-021-01546-0)

## See also

[`simulate_design()`](https://pablobernabeu.github.io/pilotr/r/reference/simulate_design.md)
to simulate from the recovered specification,
[`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md)
to run the power analysis it was read off the fit for, and
[`model_formula()`](https://pablobernabeu.github.io/pilotr/r/reference/model_formula.md)
for the analysis model a specification implies.

## Examples

``` r
# \donttest{
if (requireNamespace("lme4", quietly = TRUE)) {
  # Stand in for a pilot study: simulate a small design and fit it.
  pilot <- build_spec(list(name = "pilot", seed = 1, design_kind = "within",
    include_items = TRUE, n_subject = 30, n_item = 20, factor_name = "cond",
    lev1 = "a", lev2 = "b", intercept = 6, effect = 0.05,
    subj_int_sd = 0.12, subj_slope_sd = 0, subj_corr = 0,
    item_int_sd = 0.08, item_slope_sd = 0, item_corr = 0,
    family = "shifted_lognormal", resp_name = "", sigma = 0.3, shift = 200))
  fit <- lme4::lmer(model_formula(pilot), data = model_data(pilot, simulate_design(pilot)))

  # Read the design back off the fit, scaled up to the sample size being planned.
  spec <- spec_from_model(fit, n_subject = 60, n_item = 40)
  spec$random$subject$intercept_sd
  attr(spec, "column_kinds")
  head(simulate_design(spec))
}
#> spec_from_model() read 'effect' as categorical factor. A fitted model does not record which of its numeric columns were coded factors, so please check that reading.
#>   subject item effect_level       .y
#> 1       1    1         -0.5 5.762494
#> 2       1    1          0.5 6.435940
#> 3       1    2         -0.5 6.306309
#> 4       1    2          0.5 6.312903
#> 5       1    3         -0.5 5.963177
#> 6       1    3          0.5 6.156628
# }
```
