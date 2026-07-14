# Power and design analysis

`pilotr` turns a specification into evidence for study planning: simulate from the ground
truth, fit the analysis model and summarise across replicates. Alongside power it reports the
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
from math import sqrt
from pilotr import power_curve

n_sims = 200
curve = power_curve(spec, subject_ns=[16, 32, 48, 64, 96, 128], n_sims=n_sims)
ns = [p["n_subject"] for p in curve]
pw = [p["power"] for p in curve]
# Each power estimate is a proportion over n_sims replicates, so it carries a
# binomial Monte Carlo standard error. The shaded band is the 95% interval.
se = [sqrt(p * (1 - p) / n_sims) for p in pw]
lo = [max(0, p - 1.96 * s) for p, s in zip(pw, se)]
hi = [min(1, p + 1.96 * s) for p, s in zip(pw, se)]

fig, ax = plt.subplots(figsize=(6, 3.4))
ax.axhline(0.8, ls="--", color="grey")
ax.fill_between(ns, lo, hi, color=BLUE, alpha=0.15)
ax.plot(ns, pw, "-o", color=BLUE)
ax.set_ylim(0, 1)
ax.set_xlabel("$N$ subjects")
ax.set_ylabel("Power")
print(show(fig))
```

The shaded band shows the Monte Carlo interval, the binomial standard error of each
power estimate over the `n_sims` replicates widened to a 95% envelope.

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

`n_converged` reports how many replicates the mixed model actually fit. Convergence problems
are common in small crossed designs, so this is a useful diagnostic in its own right, and it is
the denominator of `power`, the significant proportion among the converged replicates.

Even at this tiny `n_sims`, the fixed effect is recovered (`mean_estimate` is close to
`true_effect`). One caveat: the statsmodels variance-component fit overstates random-slope
variance, so the power estimate is conservative for random-slope designs. The R package's
`lme4`-based `power_mixed` is the reference. Data generation is identical across the two
languages, and the discrepancy is in the estimator alone.

`power_mixed` runs pilotr's own simulation loop over the portable specification rather
than wrapping an existing power package. It covers territory pioneered by
[simr](https://doi.org/10.1111/2041-210x.12504) (Green and MacLeod, 2016) and
[mixedpower](https://doi.org/10.3758/s13428-021-01546-0) (Kumle, Vo and Draschkow, 2021).
pilotr differs in being driven by the portable cross-language specification, in reporting
Type S and Type M errors alongside power and in built-in parallelisation.

## Parallel execution

Every analysis on this page takes a `workers` argument that spreads the Monte Carlo
replicates across local processes with `concurrent.futures`. Each replicate seeds the
shared cross-language RNG from its own index, so the results are identical to a serial
run whatever the worker count, and parallelisation costs nothing in reproducibility. The
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
