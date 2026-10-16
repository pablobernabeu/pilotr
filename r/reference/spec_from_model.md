# Build a design specification from a fitted mixed model

Read a design specification off a linear mixed model already fitted with
[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) or
[`lmerTest::lmer()`](https://rdrr.io/pkg/lmerTest/man/lmer.html). The
fixed effects, the random-effect standard deviations and correlations,
and the residual standard deviation are taken from the fit, and the
numbers of subjects and items may be raised at the same time, so that a
pilot study or a published model becomes the starting point of a power
analysis. The fixed effects of a pilot are its point estimates, and a
study powered for them tends to be underpowered. A message names them on
every call, and `sesoi` or `effects = "safeguard"` replaces them (see
Details). Requires the `lme4` package.

## Usage

``` r
spec_from_model(
  fit,
  name = NULL,
  seed = 1,
  n_subject = NULL,
  n_item = NULL,
  family = NULL,
  round = NULL,
  effects = c("estimate", "safeguard"),
  safeguard_level = 0.6,
  sesoi = NULL,
  counterbalance = TRUE
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

- effects:

  Where the coefficients of the returned specification come from.
  `"estimate"`, the default, takes the fit's point estimates, as earlier
  versions did. `"safeguard"` takes the bound of each coefficient's
  two-sided `safeguard_level` interval nearer zero, and stops with an
  error when that interval includes zero (see Details). Either way, a
  coefficient named in `sesoi` takes the value given there, and the
  intercept is the fit's.

- safeguard_level:

  The confidence level of the intervals behind `effects = "safeguard"`,
  a single number greater than 0 and less than 1. The default of 0.6
  follows Perugini, Gallucci and Costantini (2014).

- sesoi:

  `NULL`, the default, or a named numeric vector that sets coefficients
  to a smallest effect size of interest on the scale of the fit, as in
  `sesoi = c(cond = 0.03)`. The names are the returned specification's
  coefficient keys, which write an interaction as `a:b`, and a name that
  is not among them is refused.

- counterbalance:

  Whether a factor that the pilot counterbalanced sets `counterbalance`
  in the returned specification, which then declares specification
  version 0.4 (see Details). `TRUE`, the default, plans the pilot's
  design, unless its subjects each saw only some of the items or its
  model has a term named `list`. `FALSE` returns a fully crossed
  specification that releases reading version 0.3 accept, with a warning
  that the power computed from it is too high.

## Value

A design specification as a nested list, carrying the lowest
`spec_version` that reads it, 0.4 when it sets `counterbalance` and 0.3
otherwise, and validated with
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
`"factor"` or `"predictor"`. A third, `effects`, is a list recording
where the coefficients came from. Its `method` is the value of
`effects`, its `level` the `safeguard_level` applied, or `NA` with
`"estimate"`, and its `sesoi` the keys of the coefficients set from
`sesoi`. All three are attributes rather than fields so that the
specification itself stays within the portable schema.

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

The fixed effects are another matter. They are the pilot's point
estimates, which a small pilot measures imprecisely, and a pilot
followed up because its effect looked promising has on average
overestimated that effect. A study powered for such an estimate
therefore tends to have less power than planned (Albers and Lakens,
2018). The variance components are the part of a pilot worth carrying
over. The effect is better set to a smallest effect size of interest
(SESOI), the smallest effect of theoretical or practical importance. For
an effect keyed `cond`, as in the `crossed_mixed_rt` example, pass
`sesoi = c(cond = 0.03)`, or set `spec$fixed$coefficients$cond <- 0.03`
on the returned specification. Every call reports in a message which
coefficients it took from the fit, with their values.

When no SESOI can be argued for, `effects = "safeguard"` takes safeguard
values (Perugini, Gallucci and Costantini, 2014). Each coefficient not
named in `sesoi` becomes the bound of its two-sided `safeguard_level`
interval nearer zero, by default that of the 60% interval. The interval
is the estimate plus or minus a quantile times the standard error from
`vcov(fit)`. For an `lme4` fit, the quantile is normal,
`qnorm((1 + safeguard_level) / 2)`, which is 0.84 at the default. The
normal quantile treats the standard error as known, so the interval is
too narrow when the subjects or items are few. For a `lmerTest` fit, the
quantile is a t quantile on the Satterthwaite degrees of freedom that
`lmerTest` reports, which allow for the uncertainty in the variance
components, although only approximately. When `lmerTest` cannot give
them for every coefficient, all coefficients take the normal quantile,
and the message names the distribution used. An interval that includes
zero stops the call with an error naming the effect. The pilot then
leaves even the sign of the effect open, and a value of zero would plan
a study of no effect at all. Such an effect needs a value in `sesoi`. By
default, the conditions that
[`design_conditions()`](https://pablobernabeu.github.io/pilotr/r/reference/design_conditions.md)
builds for
[`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md)
include one with every named effect at zero, which shows how often a
design declares an effect that is not there.

Whatever their source, the coefficients of the specification are the
effects its data are simulated with. A named numeric `focal` in
[`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md)
changes only the true values against which Type S and Type M are
computed. The simulated effect stays as it was, so powering a study for
another effect means changing the specification.

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

Two or more factors placed between the same unit become one factor whose
levels are their cells, and a message says so. The first-listed factor
varies slowest, and each label joins the components' labels with `"."`.
Two between factors over one unit need specification version 0.4, which
allocates them jointly. Up to version 0.3, each was assigned to blocks
of its unit on its own, which confounded their effects, and
[`validate_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/validate_spec.md)
refuses them in a specification declaring an earlier version. The factor
of cells allocates the units as version 0.4 allocates the separate
factors, and every release that reads version 0.3 reads it, so the
returned specification declares 0.3. The combined factor carries every
component's contrast columns under the names the fit gave them, so each
coefficient the fit estimated, interactions included, still names a
column. Every cell then receives units, in equal numbers whenever the
number of units is a multiple of the number of cells.

Interactions need one further step. lme4 writes the interaction of two
model-frame columns as `a:b`, which is already the specification's
convention, but
[`model_data()`](https://pablobernabeu.github.io/pilotr/r/reference/model_data.md)
gives an interaction its own product column named `a_b`, so a
specification built from a fit of pilotr's own modelling data would
otherwise acquire a spurious independent term. Such a column is
recognised by checking the product identity in the data, which both
avoids mistaking a column called `z_freq` for an interaction and settles
where to split a name with several underscores. A recognised product
column is reported and re-keyed to `a:b`, and it contributes no factor
or predictor of its own.

Grouping factors named `subject` and `item` are used as they stand.
Otherwise the one with the most levels becomes `subject`, the largest
remaining factor that is crossed with it becomes `item`, and every other
grouping factor becomes an extra `random` entry with the `over` and `n`
fields that pilotr's additional grouping factors take, `over` being
decided by which unit the factor partitions. Any renaming is reported by
a message and recorded in the `group_mapping` attribute of the result.
When each subject saw only some of the items, the item unit gains a
`per_subject` count, so that the recovered design keeps the partial
crossing of the original.

A factor that varies within both units is crossed with every
subject-item pair unless it sets `counterbalance`, a field of
specification version 0.4. A pilot run in counterbalanced lists, as
priming studies usually are, showed each subject each item at one level
only. The levels rotated across lists of subjects, so that every item
appeared at each level. When every subject-item pair of the pilot holds
one level of a factor, and every subject and every item holds two or
more, the returned specification sets `counterbalance` to `TRUE` on that
factor and declares version 0.4. Each subject then meets each item at
one level, the levels rotating over subjects and items as described
under Counterbalanced lists in the specification (`spec/SPEC.md`). The
simulated data name each subject's list in a column `list`. Factors
rotated together, as two two-level factors over four lists, are all
counterbalanced, and a message names them. Releases that read no version
above 0.3 refuse such a specification.

`counterbalance = FALSE` returns the crossed specification of earlier
versions, which those releases read, with a warning. It has twice the
pilot's observations per subject for a two-level factor, three times for
a three-level one and four times for two two-level factors rotated
together. The power computed from it is therefore too high. Take the
variance components of the `crossed_mixed_rt` example, with 30 subjects,
24 items and an effect of 0.05 on the log scale. By the formulas of
Judd, Westfall and Kenny (2017), power is then about 0.80 for the
crossed design and 0.55 for the counterbalanced one. The counterbalanced
study needs about 1.8 times the subjects to match the crossed one. A
release that reads no version above 0.3 takes a two-level design written
with a list factor between subjects and an item-set factor between
items, as in the two-list example among the worked encodings of the
specification. The crossed specification and the warning also come back,
whatever `counterbalance` says, when each subject of the pilot saw only
some of the items. The specification then keeps `per_subject`, which
version 0.4 does not combine with `counterbalance`, since that field
needs every subject to see every item. They come back as well when the
pilot's model has a column or grouping factor named `list`, as a model
that fits the list beside the condition does. Version 0.4 gives that
name to the column holding each subject's list, and no other column may
share it. The warning then asks for the term to be renamed in the
pilot's data and the model refitted.

Alternatively, keep the crossed specification and drop the rows that a
counterbalanced study would not observe, through the `prep` argument of
[`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md)
and
[`precision_design()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md).
That function receives each simulated data set and returns the modelling
data, so it ends by calling
[`model_data()`](https://pablobernabeu.github.io/pilotr/r/reference/model_data.md).
With `workers > 1` it runs in worker processes that have not attached
pilotr and cannot see the caller's variables, so it has to call
[`pilotr::model_data()`](https://pablobernabeu.github.io/pilotr/r/reference/model_data.md)
and carry the specification with it. A function made by another function
carries it once that function has forced its argument, as in this
rotation over two lists for the `crossed_mixed_rt` example, whose factor
is `condition`:

    counterbalance <- function(spec) {
      force(spec)
      function(d) {
        level <- match(d$condition, spec$factors[[1]]$levels) - 1L
        pilotr::model_data(spec, d[(d$subject %% 2L + d$item + level) %% 2L == 0L, ])
      }
    }
    power_mixed(spec, prep = counterbalance(spec), workers = 4)

A boundary-singular pilot fit is carried across as it stands, which
means a variance estimated at zero or a correlation estimated at exactly
plus or minus one. Those are faithful readings of the fit, and no kind
of defect, but they are also the sign that the random-effect structure
was richer than the pilot could support, and a power analysis resting on
them will inherit that (Bates et al., 2015). Widening the design, or
simplifying the structure before refitting, is the remedy.

## References

Albers, C. and Lakens, D. (2018). When power analyses based on pilot
data are biased: Inaccurate effect size estimators and follow-up bias.
*Journal of Experimental Social Psychology*, 74, 187-195.
[doi:10.1016/j.jesp.2017.09.004](https://doi.org/10.1016/j.jesp.2017.09.004)

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

Judd, C. M., Westfall, J. and Kenny, D. A. (2017). Experiments with more
than one random factor: Designs, analytic models, and statistical power.
*Annual Review of Psychology*, 68, 601-625.
[doi:10.1146/annurev-psych-122414-033702](https://doi.org/10.1146/annurev-psych-122414-033702)

Kumle, L., Vo, M. L.-H. and Draschkow, D. (2021). Estimating power in
(generalized) linear mixed models: An open introduction and tutorial in
R. *Behavior Research Methods*, 53, 2528-2543.
[doi:10.3758/s13428-021-01546-0](https://doi.org/10.3758/s13428-021-01546-0)

Perugini, M., Gallucci, M. and Costantini, G. (2014). Safeguard power as
a protection against imprecise power estimates. *Perspectives on
Psychological Science*, 9(3), 319-332.
[doi:10.1177/1745691614528519](https://doi.org/10.1177/1745691614528519)

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

  # The effect is the pilot's estimate. Plan for a smallest effect size of interest, here 0.03
  # on the log scale, or for the bound of the estimate's 60% interval nearer zero.
  spec_from_model(fit, n_subject = 60, n_item = 40, sesoi = c(effect = 0.03))$fixed
  spec_from_model(fit, n_subject = 60, n_item = 40, effects = "safeguard")$fixed
}
#> spec_from_model() read 'effect' as categorical factor. A fitted model does not record which of its numeric columns were coded factors, so please check that reading.
#> spec_from_model() took 'effect' = 0.0669 from the fit, as the pilot estimated it. A pilot's estimates are noisy, and a study powered for them tends to be underpowered. Set a smallest effect size of interest with `sesoi`, or use effects = "safeguard" (see ?spec_from_model).
#> spec_from_model() read 'effect' as categorical factor. A fitted model does not record which of its numeric columns were coded factors, so please check that reading.
#> spec_from_model() set 'effect' = 0.03 from `sesoi`.
#> spec_from_model() read 'effect' as categorical factor. A fitted model does not record which of its numeric columns were coded factors, so please check that reading.
#> spec_from_model() took 'effect' = 0.0526 from the fit as a safeguard value, the bound nearer zero of its 60% interval, computed with the normal distribution. The pilot estimated it at 0.0669.
#> $intercept
#> [1] 6.002566
#> 
#> $coefficients
#> $coefficients$effect
#> [1] 0.05256394
#> 
#> 
# }
```
