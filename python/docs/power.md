# Power and design analysis

`pilotr` turns a specification into evidence for study planning: simulate from the ground
truth, fit the analysis model and summarise across replicates. Alongside power it reports the
Type S (sign) and Type M (magnitude, or exaggeration) errors of Gelman and Carlin (2014),
computed over the significant replicates.

## Two-group Gaussian power

`power` handles the two-group Gaussian design with one row per subject, using a two-sample
t-test. It needs `scipy` (install the `power` extra):

```python exec="true" session="pow"
import sys; sys.path.insert(0, "docs")
from _exec import table, show, BLUE
```

```python exec="true" source="material-block" session="pow"
from pilotr import power

spec = {
    "name": "d", "seed": 1,
    "units": {"subject": {"n": 64}},
    "factors": [{"name": "group", "levels": ["a", "b"],
                 "contrasts": {"effect": [-0.5, 0.5]}, "between": "subject"}],
    "fixed": {"intercept": 100, "coefficients": {"effect": 5}},
    "response": {"family": "gaussian", "name": "score", "sigma": 10},
}

res = power(spec, n_sims=500)
print(table([{k: res[k] for k in (
    "power", "type_s", "type_m", "true_effect", "mean_estimate"
)}]))
```

At roughly 50% power the Type M ratio is well above 1. Conditional on significance the
estimated effect is exaggerated, even though the average estimate over all replicates is
unbiased. This is the statistical-significance filter that design analysis is meant to expose.

Power is a proportion over the `n_sims` replicates, so it comes with its Monte Carlo error.
`power_mcse` is its standard error, and `power_lo` and `power_hi` bound its 95% Wilson interval,
which stays between 0 and 1 and does not shrink to a point when power is 0 or 1. Type S and
Type M rest on the significant replicates alone, so they settle more slowly than power.
`type_s_lo` and `type_s_hi` bound a Wilson interval for Type S, whose plain standard error would
be zero whenever no sign is wrong, and `type_m_mcse` is the standard error of Type M.

```python exec="true" source="material-block" session="pow"
print(table([{k: res[k] for k in (
    "power_mcse", "power_lo", "power_hi", "type_s_lo", "type_s_hi", "type_m_mcse"
)}]))
```

The t-test treats every row as an independent observation, which is valid only when each
subject contributes one row and no rows share a cluster. Take 30 subjects crossed with 20
items, with by-subject and residual standard deviations of 1 and a by-item one of 0.3. With no
true effect, a row-by-row test is significant in about half of all replicates, against the
nominal 5%. `power` therefore refuses a design with an item unit, a within factor or a
grouping factor besides `subject`, and its message names what it found. The R package's
`power_mixed()` fits the model that such a design implies. Continuous predictors and
by-subject random effects are allowed, since with one row per subject they vary independently
from row to row.

## Directional predictions

A preregistration that predicts the direction of an effect can test it one-sided. In `power`,
`power_curve` and `power_mixed`, `alternative="greater"` tests whether the effect is positive
and `alternative="less"` whether it is negative. The values are spelt as in the R package, so
`"two.sided"` is the default, and scipy's spelling `"two-sided"` is refused. Here the second
group is predicted to score higher than the first.

```python exec="true" source="material-block" session="pow"
one_sided = power(spec, n_sims=500, alternative="greater")
print(table([{k: one_sided[k] for k in ("alternative", "power", "type_s", "type_m")}]))
```

The one-sided p-value comes from the same t as the two-sided one, and it is half the two-sided
p-value when the estimate lies in the predicted direction. A one-sided test at α = 0.05
therefore rejects in exactly the replicates in which the two-sided test at 0.10 rejects with the
estimate in that direction. Its power exceeds the two-sided power above. Type S and Type M
are taken over the replicates significant in the predicted direction. Type S is therefore 0 by
construction, and Type M is the exaggeration among the replicates that would support the
prediction. `power_mixed` takes its one-sided p-value from the same Wald z, and its `backend`
names the direction. Results from earlier versions were two-sided. In the direction of the true
effect, a one-sided power at α is the two-sided power at 2α times one minus the Type S at 2α.

## Power over sample size

`power_curve` sweeps the number of subjects and reports power at each:

```python exec="true" source="material-block" html="true" session="pow"
import matplotlib.pyplot as plt
from pilotr import power_curve, target_n

n_sims = 200
curve = power_curve(spec, subject_ns=[16, 32, 48, 64, 96, 128, 160, 192], n_sims=n_sims)
ns = [p["n_subject"] for p in curve]
pw = [p["power"] for p in curve]
# Each power estimate is a proportion over n_sims replicates, and each record
# carries its 95% Wilson interval, which the shaded band draws.
lo = [p["power_lo"] for p in curve]
hi = [p["power_hi"] for p in curve]

# The sample size the curve solves to, drawn as a vertical band. The next
# section explains the call.
solved = target_n(curve, target=0.8)

fig, ax = plt.subplots(figsize=(6, 3.4))
ax.axhline(0.8, ls="--", color="grey")
ax.axvspan(solved["lo"], solved["hi"], color="grey", alpha=0.12)
ax.axvline(solved["value"], ls="--", color="grey")
ax.fill_between(ns, lo, hi, color=BLUE, alpha=0.15)
ax.plot(ns, pw, "-o", color=BLUE)
ax.set_ylim(0, 1)
ax.set_xlabel("$N$ subjects")
ax.set_ylabel("Power")
print(show(fig))
```

The shaded horizontal band along the curve is the Monte Carlo interval of each power estimate,
the 95% Wilson interval over its `n_sims` replicates that `power_curve` returns as `power_lo`
and `power_hi`. The vertical band is the solved sample size, with its own interval.

Each point of the curve draws its own replicates. Point j of a sweep over J sizes takes
`replicate_seeds(seed, J)[j]` as its seed, so no two points share their random numbers, as
`target_n` assumes when it fits them. Up to pilotr 0.3.0 every point took the specification's own
seed, which `seeds="shared"` still does.

## The sample size the curve implies

The point of a power curve is the sample size at which power reaches the target, and that is
the number a preregistration quotes. Judging the crossing by eye gives a bare figure over
points whose Monte Carlo intervals overlap. `target_n` fits the curve and inverts the fit, so
the crossing arrives with the uncertainty that a simulated curve carries.

```python exec="true" source="material-block" session="pow"
solved = target_n(curve, target=0.8)
# n, n_lo and n_hi round up to whole subjects; lo and hi are the bounds before rounding.
print(table([{"target": solved["target"], "n": solved["n"],
              "n_lo": solved["n_lo"], "n_hi": solved["n_hi"],
              "lo": round(solved["lo"], 1), "hi": round(solved["hi"], 1)}]))
```

This design has a closed-form answer to check against. With an effect of 5 and a residual
standard deviation of 10, `power.t.test` in R puts 0.80 power at 127.5 subjects in total, which
falls inside the interval above. The interval spans some twenty subjects, the honest report at
200 replicates a point, and raising `n_sims` narrows it.

Nothing is extrapolated. A curve that never reaches the target within the sizes swept is
refused, and the refusal reports the range the sweep did cover. Had the sweep above stopped at
128 subjects, where power was still short of 0.80, that is what would have happened.

`solve_curve` is the general form. It takes any curve with a decision rate against a swept
value, so an effect-size sweep or a sweep over items is solved the same way, with
`transform="identity"` for an effect-size axis.

## Crossed mixed-effects power

`power_mixed` fits a crossed mixed model with `statsmodels` (needs `statsmodels` and
`pandas`), and the model is the same whatever the specification declares. In statsmodels'
notation it is `yv ~ cc`, which regresses the response on the first contrast of the design's
single within factor, with independent by-subject and by-item random intercepts and slopes on
that contrast. The response is analysed as `log(y)` for the `lognormal` family and as
`log(y - shift)` for `shifted_lognormal`, whose coefficients are on the log scale, as in the R
package. Every other family is analysed on its own scale. For `bernoulli`, `poisson`, `ordinal`
and `beta` responses, that is not the scale the coefficients are written on, which is the logit
scale, or the log scale for `poisson`. A linear model of a 0/1 response estimates a difference in
probability, for instance, and setting it against the true coefficient says nothing about
exaggeration. For these families, `power_mixed` warns and returns `mean_estimate` and `type_m`
as `nan`. `type_s` compares signs alone, which the link preserves, and power is the rate at
which the linear model's test rejects.

Only that contrast is tested, and a specification that declares anything the model leaves out
draws a warning that names each such term. In the fixed part, that is a non-zero coefficient
other than the tested contrast, a between factor or a predictor. In the random part, it is a
slope on another term, a correlation between random effects or a grouping factor besides
`subject` and `item`. The R package's `power_mixed()` fits the model the specification implies
and tests every coefficient.

The design below has one within-subject factor and crossed by-subject and by-item random
intercepts and slopes. It also declares a correlation between each unit's intercept and slope,
which the independent components leave out, so the call below warns about the correlations.

```python exec="true" source="material-block" session="pow"
from pilotr import power_mixed

spec_mixed = {
    "name": "priming", "seed": 1,
    "units": {"subject": {"n": 12}, "item": {"n": 8}},
    "factors": [{"name": "condition", "levels": ["related", "unrelated"],
                 "contrasts": {"cond": [-0.5, 0.5]},
                 "vary_within": ["subject", "item"]}],
    "fixed": {"intercept": 6, "coefficients": {"cond": 0.1}},
    "random": {
        "subject": {"intercept_sd": 0.12, "slopes": {"cond": 0.04},
                    "correlations": {"intercept~cond": 0.2}},
        "item": {"intercept_sd": 0.08, "slopes": {"cond": 0.02},
                 "correlations": {"intercept~cond": -0.1}},
    },
    "response": {"family": "shifted_lognormal", "name": "RT",
                 "sigma": 0.3, "shift": 200},
}

# A tiny replicate count keeps the docs build fast. Use 200 or more for real planning.
res = power_mixed(spec_mixed, n_sims=12)
print(table([{k: res[k] for k in (
    "n_attempted", "n_returned", "n_converged", "n_singular", "n_warning"
)}]))
print(table([{k: res[k] for k in (
    "power", "true_effect", "mean_estimate", "type_s", "type_m"
)}]))
```

The priming design above is fully crossed, so each of the 12 subjects sees all 8 items in both
conditions, 192 trials in all. A counterbalanced study, in which each subject sees each item in
one condition only, has half the observations and lower power. From specification version 0.4,
the specification writes it with `"counterbalance": true` on the condition, which rotates the
conditions over subjects and items (see
[counterbalanced lists](specification.md#counterbalanced-lists)). The condition stays a within
factor, so `power_mixed` here analyses the design as it does the crossed one. Releases that read
no version above 0.3 take the two-list encoding among the
[worked encodings](specification.md#worked-encodings), which `simulate` produces exactly as the R
package does. That encoding has no within factor, so `power_mixed` here refuses it, and its
power comes from the R package's `power_mixed()`.

The first table counts the fits as the R package does. `n_returned` counts the replicates that
returned an estimate and a p-value, and `power` is the significant proportion among them.
`n_singular` counts fits with a variance component on the boundary, judged by lme4's test of a
standard deviation below 1e-4 of the residual one. `n_warning` counts fits whose optimisers did
not converge or whose Hessian was not positive definite. A component at zero can fail
statsmodels' check of the Hessian by itself, so in a singular fit that notice is left to
`n_singular`. `n_converged` counts the fits that are neither singular nor warned. Given the
model fitted here as its `formula`, the R package's `power_mixed()` reports the same five counts
on this design. Here nearly every fit is singular, since 12 subjects and 8 items cannot support
all four variance components when two of them are as small as these slopes. Such fits stay in
`power`, as in the R package, because their fixed-effect estimates remain usable. statsmodels
also notes that a fit may be on the boundary whenever a variance component is below 0.01. On
the log scale every fit draws that notice, so it counts towards none of these.

Every replicate here is analysed with the one model, singular or not. An analysis plan may
instead preregister a simplification rule, such as dropping the slopes whenever a fit is
singular, and such a rule has power and a Type I error rate of its own. It is simulated only in
the R package, whose `power_mixed()` and `precision_design()` take a ladder of models as their
`formula` and report the share of the replicates analysed with each model.

Even at this tiny `n_sims`, the fixed effect is recovered (`mean_estimate` is close to
`true_effect`). Each replicate is fitted by REML with Powell's method, and L-BFGS takes over if
Powell's does not converge. statsmodels' default chain of optimisers stopped short of the
optimum in most crossed fits, a failure that earlier versions of this page described as
statsmodels overstating random-slope variance. Over the first 12 replicates of the crossed
reaction-time example, the standard error of the effect now matches the one `lme4` gives for the
same independent components to four decimal places in 11. The same 6 replicates are significant
in both. A replicate whose Powell fit raises is fitted again with the default chain and counts
in `n_warning`.

The p-values are Wald z tests. A Wald z treats the estimate over its standard error as normal
and so ignores the uncertainty in the variance components. With few subjects or items, it
rejects more readily than a test with Satterthwaite's degrees of freedom
([Luke, 2017](https://doi.org/10.3758/s13428-016-0809-y)), and Satterthwaite's is the test the
R package's `power_mixed()` applies. In small designs, the power here therefore tends to be
higher than R's on the same data. In the example above, 6 of the 12 replicates are significant,
against 4 in R. With 8 subjects, 6 items, random intercepts alone and no true effect, the test
here rejected the null hypothesis in 0.043 of 300 replicates (Monte Carlo standard error 0.012).
On the same data, `lme4`'s Satterthwaite test of the same model rejected in 0.030. Data
generation is identical across the two languages, so these differences lie in the model and the
test alone.

`power_mixed` carries its own simulation loop over the portable specification, with no other
power package underneath it. It covers territory pioneered by
[simr](https://doi.org/10.1111/2041-210x.12504) (Green and MacLeod, 2016) and
[mixedpower](https://doi.org/10.3758/s13428-021-01546-0) (Kumle, Vo and Draschkow, 2021).
pilotr differs in being driven by the portable cross-language specification, in reporting
Type S and Type M errors alongside power and in built-in parallelisation.

## Which model each family gets

The table sets each response family against the model fitted by this package's `power_mixed`
and by the R package's `power_mixed()` and `precision_design()`. The R package chooses the
fitter by the family, so that its estimates are on the scale the coefficients are written on. It
fits a design without random effects, such as a between-subjects design with one row per
subject, with `lm()` or `glm()`.

| Family | R `power_mixed()` | R `precision_design()` | Python `power_mixed` |
|---|---|---|---|
| `gaussian`, `exgaussian` | `lmerTest::lmer()`, Satterthwaite *t* | `lme4::lmer()`, Wald interval | MixedLM, Wald *z* |
| `lognormal`, `shifted_lognormal` | as above, on the log scale | as above, on the log scale | MixedLM on the log scale, Wald *z* |
| `bernoulli` | `lme4::glmer()`, binomial, Wald *z* | `lme4::glmer()`, Wald interval on the logit scale | linear model on the response scale, a warning, no mean estimate or Type M |
| `poisson` | `lme4::glmer()`, Poisson, Wald *z* | `lme4::glmer()`, Wald interval on the log scale | as for `bernoulli` |
| `ordinal`, `beta` | linear model on the response scale, a warning, no mean estimate or Type M, and no Type S for an interaction | linear model on the response scale, a warning, no ROPE decisions (the interval width stays) | as for `bernoulli` |

A Wald *z* with few subjects or items rejects more often than its nominal level
([Li and Redden, 2015](https://doi.org/10.1186/s12874-015-0026-x);
[Luke, 2017](https://doi.org/10.3758/s13428-016-0809-y)), so the power of a small accuracy or
count design is somewhat overstated in either package. For ordinal and Beta outcomes,
the R package's `generate_design_analysis()` writes a Bayesian design analysis on the logit scale
of their coefficients.

## Parallel execution

Every analysis on this page takes a `workers` argument that spreads the Monte Carlo
replicates across local processes with `concurrent.futures`. Each replicate takes its own
seed from `replicate_seeds()`, so the results are identical to a serial run whatever the
worker count, and parallelisation costs nothing in reproducibility. The
model fits dominate the running time, which makes the speed-up close to linear in the
number of processes. `power_curve` starts one process pool and reuses it across the whole
sweep.

```python
from pilotr import power, power_curve, power_mixed

power(spec, n_sims=2000, workers=8)          # same numbers as workers=1, sooner
power_curve(spec, subject_ns=[16, 32, 48, 64, 96, 128], n_sims=2000, workers=8)
power_mixed(spec_mixed, n_sims=500, workers=8)
```

When calling a parallel analysis from a script on Windows or macOS, put the call inside an
`if __name__ == "__main__":` block, the standard requirement for Python's spawn-based
multiprocessing. Interactive sessions and notebooks need no guard.

This design answers a serial bottleneck familiar from `simr::powerCurve()` in R, which
pilotr's maintainer previously worked around by splitting the sample-size grid across
separate jobs by hand and recombining the results afterwards
([Bernabeu, 2021](https://pablobernabeu.github.io/2021/parallelizing-simr-powercurve/)).
In pilotr the same gain takes one argument.
