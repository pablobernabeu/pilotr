# Response families

The same linear predictor can be mapped through any of seven response families. The family and
its parameters live in the spec's `response` block, and the response column is named by
`response["name"]`.

| Outcome | `family` | Link | Key parameters |
|---|---|---|---|
| Continuous | `gaussian` | identity | `sigma` |
| Reaction time | `shifted_lognormal` | log | `sigma`, `shift` |
| Positive continuous (e.g. reading time) | `lognormal` | log | `sigma` |
| Accuracy (0/1) | `bernoulli` | logit | none |
| Counts | `poisson` | log | none |
| Likert / ordered | `ordinal` | cumulative logit | `thresholds` |
| Proportions in (0, 1) | `beta` | logit (mean) | `phi` (precision) |

The fixed `effect` is always expressed on the scale of the linear predictor (the link scale):
the identity scale for Gaussian, the log scale for reaction times and counts, and the logit
scale for accuracy, ordinal and Beta.

## What the families look like

```python exec="true" session="fam"
import sys; sys.path.insert(0, "docs")
from _exec import table, show, BLUE, RED, GREEN
```

A large two-group draw makes each family's shape clear: Gaussian is symmetric, the shifted
lognormal has the right skew of reaction times, the plain lognormal is the same right-skewed
shape without the shift and suits per-word reading times, and the Beta is bounded in (0, 1).

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

fig, ax = plt.subplots(1, 4, figsize=(12, 3))
ax[0].hist(draw("gaussian", 100, 5, "score", sigma=10), bins=40, color=BLUE, edgecolor="none")
ax[0].set_title("Gaussian"); ax[0].set_xlabel("score")
ax[1].hist(draw("shifted_lognormal", 6, 0.1, "RT", sigma=0.3, shift=200),
           bins=40, color=RED, edgecolor="none")
ax[1].set_title("Shifted lognormal (RT)"); ax[1].set_xlabel("RT (ms)")
ax[2].hist(draw("lognormal", -1.2, 0.1, "reading_time", sigma=0.25),
           bins=40, color=RED, edgecolor="none")
ax[2].set_title("Lognormal (reading time)"); ax[2].set_xlabel("reading time")
ax[3].hist(draw("beta", 0, 0.8, "proportion", phi=8), bins=40, color=GREEN, edgecolor="none")
ax[3].set_title("Beta"); ax[3].set_xlabel("proportion")
for a in ax:
    a.set_ylabel("count")
print(show(fig))
```

## Discrete and bounded outcomes

The same machinery covers binary, count and proportion outcomes. The fixed effect moves the two
group means apart on each family's own scale:

```python exec="true" source="material-block" session="fam"
from statistics import mean
from pilotr import simulate

def group_means(family, intercept, effect, name, **resp):
    spec = {
        "name": family, "seed": 1, "units": {"subject": {"n": 4000}},
        "factors": [{"name": "group", "levels": ["control", "treatment"],
                     "contrasts": {"effect": [-0.5, 0.5]}, "between": "subject"}],
        "fixed": {"intercept": intercept, "coefficients": {"effect": effect}},
        "response": {"family": family, "name": name, **resp},
    }
    d = simulate(spec)
    return {g: mean(r[name] for r in d.rows if r["group"] == g) for g in ("control", "treatment")}

rows = []
for label, fam, ic, ef, nm, kw in [
    ("Bernoulli: P(correct)", "bernoulli", 0.0, 0.5, "accuracy", {}),
    ("Poisson: mean count", "poisson", 1.5, 0.3, "count", {}),
    ("Beta: mean proportion", "beta", 0.0, 0.8, "p", {"phi": 8}),
]:
    m = group_means(fam, ic, ef, nm, **kw)
    rows.append({"outcome": label, "control": m["control"], "treatment": m["treatment"]})
print(table(rows))
```

An ordinal (Likert) design sets the category thresholds directly; the effect shifts mass across
the categories. The category proportions by group:

```python exec="true" source="material-block" session="fam"
spec = {
    "name": "likert", "seed": 1, "units": {"subject": {"n": 4000}},
    "factors": [{"name": "group", "levels": ["control", "treatment"],
                 "contrasts": {"effect": [-0.5, 0.5]}, "between": "subject"}],
    "fixed": {"intercept": 0.0, "coefficients": {"effect": 0.8}},
    "response": {"family": "ordinal", "name": "rating", "thresholds": [-2, -0.6, 0.6, 2]},
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

Every family ships as a ready-to-run [worked example](examples.md), and the full format,
including the thresholds and the Beta precision `phi`, is on the
[Specification](specification.md) page.
