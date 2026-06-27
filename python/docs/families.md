# Response families

The same linear predictor can be mapped through any of six response families. The family and
its parameters live in the spec's `response` block, and the response column is named by
`response["name"]`.

| Outcome | `family` | Link | Key parameters |
|---|---|---|---|
| Continuous | `gaussian` | identity | `sigma` |
| Reaction time | `shifted_lognormal` | log | `sigma`, `shift` |
| Accuracy (0/1) | `bernoulli` | logit | none |
| Counts | `poisson` | log | none |
| Likert / ordered | `ordinal` | cumulative logit | `thresholds` |
| Proportions in (0, 1) | `beta` | logit (mean) | `phi` (precision) |

The fixed `effect` is always expressed on the response scale: the identity scale for Gaussian,
the log scale for reaction times and counts, and the logit scale for accuracy, ordinal, and
Beta.

## What the families look like

```python exec="true" session="fam"
import sys; sys.path.insert(0, "docs")
from _exec import table, show, BLUE, RED, GREEN
```

A large two-group draw makes each family's shape clear: Gaussian is symmetric, the shifted
lognormal has the right skew of reaction times, and the Beta is bounded in (0, 1).

```python exec="true" source="material-block" html="true" session="fam"
import matplotlib.pyplot as plt
from pilotr import simulate

def draw(family, intercept, effect, name, **resp):
    spec = {
        "name": family, "seed": 1, "units": {"subject": {"n": 4000}},
        "factors": [{"name": "group", "levels": ["control", "treatment"],
                     "contrasts": {"effect": [-0.5, 0.5]}, "between": "subject"}],
        "fixed": {"intercept": intercept, "coefficients": {"effect": effect}},
        "response": {"family": family, "name": name, **resp},
    }
    return simulate(spec).column(name)

fig, ax = plt.subplots(1, 3, figsize=(9, 3))
ax[0].hist(draw("gaussian", 100, 5, "score", sigma=10), bins=40, color=BLUE, edgecolor="white")
ax[0].set_title("Gaussian"); ax[0].set_xlabel("score")
ax[1].hist(draw("shifted_lognormal", 6, 0.1, "RT", sigma=0.3, shift=200),
           bins=40, color=RED, edgecolor="white")
ax[1].set_title("Shifted lognormal (RT)"); ax[1].set_xlabel("RT (ms)")
ax[2].hist(draw("beta", 0, 0.8, "proportion", phi=8), bins=40, color=GREEN, edgecolor="white")
ax[2].set_title("Beta"); ax[2].set_xlabel("proportion")
for a in ax:
    a.set_ylabel("count")
print(show(fig))
```

A Bernoulli accuracy design:

```python
from pilotr import simulate

spec = {
    "name": "acc", "seed": 1,
    "units": {"subject": {"n": 80}},
    "factors": [{"name": "group", "levels": ["a", "b"],
                 "contrasts": {"effect": [-0.5, 0.5]}, "between": "subject"}],
    "fixed": {"intercept": 0.0, "coefficients": {"effect": 0.5}},
    "response": {"family": "bernoulli", "name": "accuracy"},
}
set(simulate(spec).column("accuracy"))   # {0, 1}
```

An ordinal (5-point Likert) design needs four `thresholds`; a Beta design needs a precision
`phi`. Only the `response` block changes:

```python
"response": {"family": "ordinal", "name": "rating", "thresholds": [-2, -0.6, 0.6, 2]}
"response": {"family": "beta", "name": "proportion", "phi": 8}
```

The worked example specs in
[`spec/examples/`](https://github.com/pablobernabeu/pilotr/tree/main/spec/examples) cover one
design per family.
