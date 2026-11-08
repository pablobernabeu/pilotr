# Response families

The same linear predictor can be mapped through any of eight response families. The family and
its parameters are set in the specification's `response` block, and the response column is named
by `response["name"]`.

| Outcome | `family` | Link | Key parameters |
|---|---|---|---|
| Continuous | `gaussian` | identity | `sigma` |
| Reaction time | `shifted_lognormal` | log | `sigma`, `shift` |
| Positive continuous (e.g. reading time) | `lognormal` | log | `sigma` |
| Reaction time, on the response scale | `exgaussian` | identity | `sigma`, `beta` |
| Accuracy (0/1) | `bernoulli` | logit | none |
| Counts | `poisson` | log | none |
| Likert / ordered | `ordinal` | cumulative logit | `thresholds` |
| Proportions in (0, 1) | `beta` | logit (mean) | `phi` (precision) |

The fixed `effect` is always expressed on the scale of the linear predictor, the link scale.
That is the identity scale for the Gaussian and ex-Gaussian families, the log scale for
lognormal reaction times, reading times and counts, and the logit scale for accuracy, ordinal
and proportion outcomes.

The ex-Gaussian is a normal plus an exponential, mean-centred by subtracting the mean of the
exponential so that the linear predictor remains the mean of the response. That is brms's
`exgaussian(mu, sigma, beta)` parameterisation, so a specification and the model fitted to it
agree on what the intercept means. A shifted lognormal is not a substitute, because
`model_data` logs the response back and leaves a symmetric residual on the analysis scale. A
specification that uses the ex-Gaussian must declare a `spec_version` of `"0.3"` or later, since
a 0.2 implementation would read it differently.

Every simulated Beta response lies strictly inside (0, 1), as a Beta regression requires. Near a
ceiling or a floor, a draw can round to exactly 1 or 0 in double precision, and it then becomes
the nearest double inside the interval. When a Beta specification sets `round`, the setting must
be at least 1. A response rounded onto 0 or 1 then becomes the nearest value inside the interval
at that many decimals, 0.001 or 0.999 at three. A linear predictor beyond about 36.7, or below
about −745, leaves a shape parameter of zero, and a lognormal response whose
`exp(eta + sigma * z)` overflows has no value to return. `simulate` raises a `ValueError` in both
cases, with the same message as the R package.

## What the families look like

```python exec="true" session="fam"
import sys; sys.path.insert(0, "docs")
from _exec import table, show, BLUE, RED, GREEN
```

A large two-group draw shows the shape of each family. The Gaussian is symmetric. The shifted
lognormal has the right skew of reaction times, and the plain lognormal, which suits per-word
reading times, has the same shape without the shift. The Beta is bounded in (0, 1).

```python exec="true" source="material-block" html="true" session="fam"
import matplotlib.pyplot as plt
from pilotr import simulate

def draw(family, intercept, effect, name, **resp):
    spec = {
        "name": family, "seed": 1, "units": {"subject": {"n": 4000}},
        "factors": [{"name": "group", "levels": ["control", "treatment"],
                     "contrasts": {"effect": [-0.5, 0.5]},
                     "between": "subject"}],
        "fixed": {"intercept": intercept, "coefficients": {"effect": effect}},
        "response": {"family": family, "name": name, **resp},
    }
    return simulate(spec).column(name)

fig, ax = plt.subplots(1, 4, figsize=(12, 3))
ax[0].hist(
    draw("gaussian", 100, 5, "score", sigma=10),
    bins=40, color=BLUE, edgecolor="none")
ax[0].set_title("Gaussian"); ax[0].set_xlabel("score")
ax[1].hist(draw("shifted_lognormal", 6, 0.1, "RT", sigma=0.3, shift=200),
           bins=40, color=RED, edgecolor="none")
ax[1].set_title("Shifted lognormal (RT)"); ax[1].set_xlabel("RT (ms)")
ax[2].hist(draw("lognormal", -1.2, 0.1, "reading_time", sigma=0.25),
           bins=40, color=RED, edgecolor="none")
ax[2].set_title("Lognormal (reading time)"); ax[2].set_xlabel("reading time")
ax[3].hist(
    draw("beta", 0, 0.8, "proportion", phi=8),
    bins=40, color=GREEN, edgecolor="none")
ax[3].set_title("Beta"); ax[3].set_xlabel("proportion")
for a in ax:
    a.set_ylabel("count")
print(show(fig))
```

## Discrete and bounded outcomes

Binary, count and proportion outcomes are simulated in the same way. The fixed effect moves the
two group means apart on each family's own scale:

```python exec="true" source="material-block" session="fam"
from statistics import mean
from pilotr import simulate

def group_means(family, intercept, effect, name, **resp):
    spec = {
        "name": family, "seed": 1, "units": {"subject": {"n": 4000}},
        "factors": [{"name": "group", "levels": ["control", "treatment"],
                     "contrasts": {"effect": [-0.5, 0.5]},
                     "between": "subject"}],
        "fixed": {"intercept": intercept, "coefficients": {"effect": effect}},
        "response": {"family": family, "name": name, **resp},
    }
    d = simulate(spec)
    return {
        g: mean(r[name] for r in d.rows if r["group"] == g)
        for g in ("control", "treatment")
    }

rows = []
for label, fam, ic, ef, nm, kw in [
    ("Bernoulli: P(correct)", "bernoulli", 0.0, 0.5, "accuracy", {}),
    ("Poisson: mean count", "poisson", 1.5, 0.3, "count", {}),
    ("Beta: mean proportion", "beta", 0.0, 0.8, "p", {"phi": 8}),
]:
    m = group_means(fam, ic, ef, nm, **kw)
    rows.append({"outcome": label, "control": m["control"],
                 "treatment": m["treatment"]})
print(table(rows))
```

An ordinal (Likert) design sets the category thresholds directly. The effect shifts mass across
the categories, as the category proportions by group show.

```python exec="true" source="material-block" session="fam"
spec = {
    "name": "likert", "seed": 1, "units": {"subject": {"n": 4000}},
    "factors": [{"name": "group", "levels": ["control", "treatment"],
                 "contrasts": {"effect": [-0.5, 0.5]}, "between": "subject"}],
    "fixed": {"intercept": 0.0, "coefficients": {"effect": 0.8}},
    "response": {"family": "ordinal", "name": "rating",
                 "thresholds": [-2, -0.6, 0.6, 2]},
}
d = simulate(spec)
cats = sorted(set(d.column("rating")))
rows = []
for g in ("control", "treatment"):
    vals = [r["rating"] for r in d.rows if r["group"] == g]
    row = {"group": g}
    row.update({f"P({c})": vals.count(c) / len(vals) for c in cats})
    rows.append(row)
print(table(rows))
```

Six of the eight families, all but `bernoulli` and `exgaussian`, appear in a ready-to-run
[worked example](examples.md). The [Specification](specification.md) page defines the format,
including the thresholds and the Beta precision `phi`.
