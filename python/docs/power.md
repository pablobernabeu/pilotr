# Power and design analysis

`pilotr` turns a specification into evidence for study planning: simulate from the ground
truth, fit the analysis model, and summarize across replicates. Alongside power it reports the
Type S (sign) and Type M (magnitude, or exaggeration) errors of Gelman and Carlin (2014),
computed over the significant replicates.

## Two-group Gaussian power

`power` handles the two-group Gaussian design with a two-sample t-test. It needs `scipy`
(install the `power` extra):

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
print(table([{k: res[k] for k in ("power", "type_s", "type_m", "true_effect", "mean_estimate")}]))
```

At roughly 50% power the Type M ratio is well above 1. Conditional on significance the
estimated effect is exaggerated, even though the average estimate over all replicates is
unbiased. This is the statistical-significance filter that design analysis is meant to expose.

## Power over sample size

`power_curve` sweeps the number of subjects so you can read off where power crosses a target:

```python exec="true" source="material-block" html="true" session="pow"
import matplotlib.pyplot as plt
from pilotr import power_curve

curve = power_curve(spec, subject_ns=[16, 32, 48, 64, 96, 128], n_sims=200)
ns = [p["n_subject"] for p in curve]
pw = [p["power"] for p in curve]

fig, ax = plt.subplots(figsize=(6, 3.4))
ax.axhline(0.8, ls="--", color="grey")
ax.plot(ns, pw, "-o", color=BLUE)
ax.set_ylim(0, 1)
ax.set_xlabel("$N$ subjects")
ax.set_ylabel("Power")
print(show(fig))
```

## Crossed mixed-effects power

`power_mixed` fits a crossed mixed model with `statsmodels` (needs `statsmodels` and
`pandas`). The design has one within-subject factor and crossed by-subject and by-item random
intercepts and slopes:

```python exec="true" source="material-block" session="pow"
from pilotr import power_mixed

spec_mixed = {
    "name": "priming", "seed": 1,
    "units": {"subject": {"n": 12}, "item": {"n": 8}},
    "factors": [{"name": "condition", "levels": ["related", "unrelated"],
                 "contrasts": {"cond": [-0.5, 0.5]}, "vary_within": "subject"}],
    "fixed": {"intercept": 6, "coefficients": {"cond": 0.1}},
    "random": {
        "subject": {"intercept_sd": 0.12, "slopes": {"cond": 0.04},
                    "correlations": {"intercept~cond": 0.2}},
        "item": {"intercept_sd": 0.08, "slopes": {"cond": 0.02},
                 "correlations": {"intercept~cond": -0.1}},
    },
    "response": {"family": "shifted_lognormal", "name": "RT", "sigma": 0.3, "shift": 200},
}

res = power_mixed(spec_mixed, n_sims=12)   # tiny for the docs; use >= 200 in practice
print(table([{k: res[k] for k in
              ("power", "n_converged", "true_effect", "mean_estimate", "type_s", "type_m")}]))
```

Even at this tiny `n_sims`, the fixed effect is recovered (`mean_estimate` is close to
`true_effect`). One caveat: the statsmodels variance-component fit overstates random-slope
variance, so the power estimate is conservative for random-slope designs. The R package's
`lme4`-based `power_mixed` is the reference; data generation is identical across the two
languages, and the discrepancy is in the estimator alone.
