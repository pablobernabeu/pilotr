# pilotr (Python)

`pilotr` simulates experimental and behavioral data from a portable JSON **design
specification**. The same specification drives the R package, a no-code web app, and this
Python package; given the same spec and seed, all three produce bit-identical data.

This is the Python documentation. For the R reference and the live no-code demo, see the
[project documentation](https://pablobernabeu.github.io/pilotr/).

## Install

The package is not yet on PyPI. Install it from a checkout of the repository:

```bash
git clone https://github.com/pablobernabeu/pilotr.git
pip install ./pilotr/python          # add the [power] extra for the scipy-based power demo
```

The generative core is pure Python with no required dependencies. `scipy` (for `power`) and
`statsmodels` with `pandas` (for `power_mixed`) are optional extras.

## Quick start

A design is a plain dictionary, or a JSON file, describing the units, factors, fixed effects,
optional random effects, and a response family. Simulate from it with `simulate`:

```python
from pilotr import simulate

spec = {
    "name": "two_group", "seed": 2024,
    "units": {"subject": {"n": 64}},
    "factors": [{"name": "group", "levels": ["control", "treatment"],
                 "contrasts": {"effect": [-0.5, 0.5]}, "between": "subject"}],
    "fixed": {"intercept": 100, "coefficients": {"effect": 5}},
    "response": {"family": "gaussian", "name": "score", "sigma": 10},
}

data = simulate(spec)
len(data)                   # 64
data.head(3)                # first rows as a list of dicts
data.to_csv("data.csv")     # write to CSV
```

To run a spec authored elsewhere (for example, one downloaded from the no-code app), load it
with `load_spec` and simulate:

```python
from pilotr import load_spec, simulate
data = simulate(load_spec("design.json"))
```

## Where to go next

- [Cross-language reproducibility](cross-language.md): how the same spec gives identical data
  in R and Python.
- [Power and design analysis](power.md): `power`, `power_curve`, and `power_mixed`.
- [Response families](families.md): Gaussian, reaction times, accuracy, counts, ordinal, Beta.
- [API reference](api.md): every public function and class.

The specification format itself is documented in
[`spec/SPEC.md`](https://github.com/pablobernabeu/pilotr/blob/main/spec/SPEC.md).
