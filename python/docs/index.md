# pilotr (Python)

`pilotr` simulates experimental and behavioural data from a portable JSON design
specification. The same specification drives the
[R package](https://pablobernabeu.github.io/pilotr/), a no-code web app and this Python package.
Given the same spec and seed, all three produce bit-identical data.

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
optional random effects and a response family. Simulate from it with `simulate`:

!!! note "These examples run live"
    The code blocks on this page (and throughout these docs) are executed when the site is
    built, so the tables and plots are real `pilotr` output. `table(...)` and `show(...)` are
    small helpers that render a result as a Markdown table or an inline figure.

```python exec="true" session="quick"
import sys; sys.path.insert(0, "docs")
from _exec import table, show, BLUE, RED, GREEN
```

```python exec="true" source="material-block" session="quick"
from pilotr import simulate

spec = {
    "name": "two_group", "seed": 2024,
    "units": {"subject": {"n": 64}},
    "factors": [{"name": "group", "levels": ["control", "treatment"],
                 "contrasts": {"effect": [-0.5, 0.5]}, "between": "subject"}],
    "fixed": {"intercept": 100, "coefficients": {"effect": 5}},
    "response": {"family": "gaussian", "name": "score", "sigma": 10},
}

data = simulate(spec)               # 64 rows
print(table(data.head(5)))          # the first rows, as a table
```

`len(data)` is the number of observations, `data.head(n)` returns the first rows as a list of
dicts, and `data.to_csv("data.csv")` writes the table to disk.

To run a spec authored elsewhere (for example, one downloaded from the no-code app), load it
with `load_spec` and simulate. Here we load one of the worked examples that ship with pilotr:

```python exec="true" source="material-block" session="quick"
from pilotr import load_spec, simulate

data = simulate(load_spec("../spec/examples/poisson_counts_between.json"))
print(table(data.head(5)))
```

## Where to go next

- [Cross-language reproducibility](cross-language.md): how the same spec gives identical data
  in R and Python.
- [Power and design analysis](power.md): `power`, `power_curve`, and `power_mixed`.
- [Response families](families.md): Gaussian, reaction times, accuracy, counts, ordinal, Beta.
- [Worked examples](examples.md): one ready-to-run design per family.
- [Specification format](specification.md): the JSON spec and the cross-language RNG contract.
- [API reference](api.md): every public function and class.
