# Get started

This guide installs pilotr and simulates a first data set from a design specification. The
specification is what everything else reads. The same JSON file drives the
[R package](https://pablobernabeu.github.io/pilotr/), the
[no-code app](https://pablobernabeu.github.io/pilotr/app/) and this Python package, and given
the same file and seed all three produce bit-identical data.

## Install

```bash
pip install pilotr             # core engine (pure Python, dependency-free)
pip install "pilotr[power]"    # + scipy, for the simulation-based power demo
pip install "pilotr[mixed]"    # + scipy, statsmodels and pandas, for crossed mixed-effects power
```

Requires Python 3.9 or later. The generative core has no dependencies. `scipy` (for `power`)
and `statsmodels` with `pandas` (for `power_mixed`) are optional extras, so an install that
only simulates stays dependency-free. For development, install from a checkout of the
repository instead:

```bash
git clone https://github.com/pablobernabeu/pilotr.git
pip install ./pilotr/python
```

## A first simulation

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

The specification is checked before anything is drawn, because several ways of getting one
wrong produce plausible data rather than an error. A mistyped coefficient key resolves to no
column and so silently sets that effect to zero, which generates exactly the data of a null
design. Pass `validate=False` to skip the check in a loop that has already validated the
specification once.

## Running a specification authored elsewhere

To run a spec authored elsewhere, for example one downloaded from the no-code app, load it with
`load_spec` and simulate. Here we load one of the worked examples that ship with pilotr, using
`pilotr_example` to find it inside the installed package:

```python exec="true" source="material-block" session="quick"
from pilotr import load_spec, pilotr_example, simulate

data = simulate(load_spec(pilotr_example("poisson_counts_between")))
print(table(data.head(5)))
```

`pilotr_example()` with no argument lists every bundled specification, one per design family.

## Where to go next

From here, [Response families](families.md) shows what each of the eight families generates and
which scale its effect is written on, [Worked examples](examples.md) gives one ready-to-run
design per family, and [Power and design analysis](power.md) turns a specification into a power
estimate with its Type S and Type M errors. [Cross-language reproducibility](cross-language.md)
explains why an R run and a Python run agree to the last bit,
[Specification format](specification.md) documents the JSON itself, and the
[API reference](api.md) lists every public function and class.
