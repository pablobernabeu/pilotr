# Cross-language reproducibility

The point of the shared specification is that one design produces the same data in R and in
Python. The two implementations never call each other; they share a numerical contract.

## The shared generator

Native random-number generators differ across ecosystems (R uses Mersenne-Twister with
inversion, NumPy uses PCG64), so a naive port would diverge. `pilotr` instead ships one
generator implemented identically in both languages.

- **Uniforms** come from L'Ecuyer's (1988) combined linear congruential generator. Every
  intermediate product stays below `2**53`, so the arithmetic is exact in IEEE-754 doubles
  and in Python integers alike.
- **Normals** use Wichura's (1988) Algorithm AS 241 inverse-CDF, the same routine R's
  `qnorm` uses, so deviates agree to full double precision.
- **Everything else** (Cholesky-correlated random effects, inverse-CDF Poisson and ordinal
  draws, Marsaglia and Tsang gamma draws for the Beta family) derives from those two through
  a documented, identical draw order.

The draw order is specified in
[`spec/SPEC.md`](https://github.com/pablobernabeu/pilotr/blob/main/spec/SPEC.md).

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

The numeric columns match to the last digit. The repository's `python/examples/parity_check.py`
runs this comparison across the worked example designs and reports the maximum difference,
which is zero up to floating-point accumulation.

A simulated crossed design has one row per subject-by-item observation. The first rows of the
worked crossed reaction-time example look like this (the same rows the R package produces):

```python exec="true" session="xl"
import sys; sys.path.insert(0, "docs")
from _exec import table, show
```

```python exec="true" source="material-block" session="xl"
from pilotr import simulate, load_spec

data = simulate(load_spec("../spec/examples/crossed_mixed_rt.json"))
print(table(data.head(6)))
```
