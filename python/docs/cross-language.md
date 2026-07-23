# Cross-language reproducibility

The point of the shared specification is that one design produces the same data in R and in
Python. The two implementations never call each other; they share a numerical contract.

## The shared generator

Native random-number generators differ across ecosystems (R uses Mersenne-Twister with
inversion, NumPy uses PCG64), so a naive port would diverge. `pilotr` instead ships one
generator implemented identically in both languages.

Uniforms come from L'Ecuyer's (1988) combined linear congruential generator. Every
intermediate product stays below `2**53`, so the arithmetic is exact in IEEE-754 doubles
and in Python integers alike. Normals use Wichura's (1988) Algorithm AS 241 inverse-CDF, the
same routine R's `qnorm` uses, so deviates agree to full double precision. Everything else
(Cholesky-correlated random effects, inverse-CDF Poisson and ordinal draws, Marsaglia and
Tsang gamma draws for the Beta family) derives from those two through a documented, identical
draw order.

The draw order is specified on the [Specification](specification.md) page.

## Determinism in Python

The same specification and seed always give the same data.

```python exec="true" source="material-block" session="xl"
from pilotr import simulate

spec = {
    "name": "demo", "seed": 2024,
    "units": {"subject": {"n": 200}},
    "factors": [{"name": "group", "levels": ["a", "b"],
                 "contrasts": {"effect": [-0.5, 0.5]}, "between": "subject"}],
    "fixed": {"intercept": 100, "coefficients": {"effect": 5}},
    "response": {"family": "gaussian", "name": "score", "sigma": 10},
}

print(simulate(spec).rows == simulate(spec).rows)
```

Changing the seed changes the draws, while the structure and the ground-truth parameters stay
fixed.

```python exec="true" source="material-block" session="xl"
spec2 = dict(spec, seed=spec["seed"] + 1)
print(simulate(spec).column("score") == simulate(spec2).column("score"))
```

## Verifying parity

Simulate the same spec and seed in both languages and compare. In Python:

```python
from pilotr import simulate, load_spec
simulate(load_spec("design.json")).to_csv("py.csv")
```

In R:

```r
library(pilotr)
write.csv(simulate_design(load_spec("design.json")), "r.csv", row.names = FALSE)
```

The two CSVs compare exactly, up to CSV number formatting. The repository's [`python/examples/parity_check.py`](https://github.com/pablobernabeu/pilotr/blob/main/python/examples/parity_check.py)
runs this comparison across the worked example designs and reports the maximum difference:
zero for every design with rounded responses, and below 1e-14 for the unrounded continuous
design, where R's CSV writer prints 15 significant digits. The script's small tolerance
(1e-6) absorbs that number formatting only, never a divergence in the generators.

A simulated crossed design has one row per subject-by-item observation. The first rows of the
worked crossed reaction-time example look like this (the same rows the R package produces):

```python exec="true" session="xl"
import sys; sys.path.insert(0, "docs")
from _exec import table, show
```

```python exec="true" source="material-block" session="xl"
from pilotr import simulate, load_spec, pilotr_example

data = simulate(load_spec(pilotr_example("crossed_mixed_rt")))
print(table(data.head(6)))
```
