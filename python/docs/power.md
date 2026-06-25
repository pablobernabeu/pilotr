# Power and design analysis

`pilotr` turns a specification into evidence for study planning: simulate from the ground
truth, fit the analysis model, and summarize across replicates. Alongside power it reports the
Type S (sign) and Type M (magnitude, or exaggeration) errors of Gelman and Carlin (2014),
computed over the significant replicates.

## Two-group Gaussian power

`power` handles the two-group Gaussian design with a two-sample t-test. It needs `scipy`
(install the `power` extra):

```python
from pilotr import power

spec = {
    "name": "d", "seed": 1,
    "units": {"subject": {"n": 64}},
    "factors": [{"name": "group", "levels": ["a", "b"],
                 "contrasts": {"effect": [-0.5, 0.5]}, "between": "subject"}],
    "fixed": {"intercept": 100, "coefficients": {"effect": 5}},
    "response": {"family": "gaussian", "name": "score", "sigma": 10},
}

power(spec, n_sims=500)
# {'power': ~0.50, 'type_s': 0.0, 'type_m': ~1.3, 'true_effect': 5.0, ...}
```

At roughly 50% power the Type M ratio is well above 1. Conditional on significance the
estimated effect is exaggerated, even though the average estimate over all replicates is
unbiased. This is the statistical-significance filter that design analysis is meant to expose.

## Power over sample size

`power_curve` sweeps the number of subjects so you can read off where power crosses a target:

```python
from pilotr import power_curve
power_curve(spec, subject_ns=[32, 48, 64, 96, 128], n_sims=200)
```

## Crossed mixed-effects power

`power_mixed` fits a crossed mixed model with `statsmodels` (needs `statsmodels` and
`pandas`). One caveat: the statsmodels variance-component fit overstates random-slope variance,
so the backend is conservative for random-slope designs. For those designs the R package's
`lme4`-based `power_mixed` is the reference; data generation is identical across the two
languages, and the discrepancy is in the estimator alone.

```python
from pilotr import power_mixed
power_mixed(spec_with_one_within_factor, n_sims=200)
```
