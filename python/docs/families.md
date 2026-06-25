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
