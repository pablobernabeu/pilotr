# pilotr (Python)

Python implementation of **pilotr** — simulate experimental and behavioral data from a
portable JSON design specification, with **bit-identical** output to the R package of the
same name. See the [project README](https://github.com/pablobernabeu/pilotr) for the full
toolkit (no-code app, R package, the specification format, and the manuscript).

## Install

```bash
pip install .            # core engine (pure Python, dependency-free)
pip install ".[power]"   # + scipy, for the simulation-based power demo
```

## Quick start

```python
from pilotr import simulate, load_spec

data = simulate("../spec/examples/between_2group_gaussian.json")
```

Given the same specification and seed, this produces the same data as the R package to full
floating-point precision.
