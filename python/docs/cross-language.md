# Cross-language reproducibility

The point of the shared specification is that one design produces the same data in R and in
Python. The two implementations never call each other. What they share is a numerical contract.

## The shared generator

Native random-number generators differ across ecosystems (R uses Mersenne-Twister with
inversion, NumPy uses PCG64), so a naive port would diverge. `pilotr` instead ships one
generator implemented identically in both languages.

Uniforms come from L'Ecuyer's (1988) combined linear congruential generator. Every
intermediate product stays below `2**53`, so the arithmetic is exact in IEEE-754 doubles
and in Python integers alike. That holds for seeds within ±(2^53 − 1), the range `validate_spec`
accepts: R reads a JSON number as a double, so beyond that range the two readers would hold
different seeds. Normals use Wichura's (1988) Algorithm AS 241 inverse-CDF, the
same routine R's `qnorm` uses. Everything else (Cholesky-correlated random effects,
inverse-CDF Poisson and ordinal draws, Marsaglia and Tsang gamma draws for the Beta family)
derives from those two through a documented, identical draw order.

The draw order is specified on the [Specification](specification.md) page.

## Where the data can differ

Given the same specification and seed, the two implementations produce identical data wherever
they share a maths library, as on Linux. Continuous integration checks every case of the parity
harness to the last bit there. IEEE-754 fixes the rounding of arithmetic but leaves that of
`exp()`, `log()` and `pow()` to the library, and on Windows Python and R do not share one.
CPython takes them from Microsoft's Universal C Runtime, and R from the MinGW-w64 toolchain it is
built with. There, a value that passes through those functions can differ in its last bits, and
every family has such a path. The figures below were measured on Windows 11 with Python 3.13 and
R 4.6.1.

- Normal draws pass through `log()` in AS 241's tail branch, which about 15% of draws take.
  About 1.4e-5 of normal draws differed, by 1 to 3 units in the last place (ulp), so Gaussian
  data can differ too.
- Responses of the lognormal families and the ex-Gaussian differed by 1 ulp.
- An unrounded Beta response differed by up to 35 ulp.
- A discrete outcome can in principle move by a whole category, although none did in 50,000
  rows of each discrete family.

One Beta case can differ on Linux too. From R 4.6.0, R's `^` computes a cube or a fourth power
by repeated multiplication, where Python calls `pow()`. A Beta response whose Gamma shape is
exactly 1/3 or 1/4 can therefore differ in its last bits on any platform. Such a shape needs a
cell whose rows share one linear predictor, as in a Beta design with no random effects, a linear
predictor of 0 and `phi` set to 0.5.

Rounding has removed every continuous difference observed so far, but it guarantees nothing,
since a last-bit difference can straddle a rounding boundary, and `round` is refused for the
discrete families. Agreement to the last bit in every family has one route: run both languages
on a platform where they share a maths library, outside that Beta case. The section Scope of the
guarantee on the [Specification](specification.md) page gives the measurements in full.

## Determinism in Python

On one platform, the same specification and seed always give the same data.

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

Across operating systems, the same holds on the terms set out above, since CPython takes
`exp()`, `log()` and `pow()` from the platform's C library. On Windows, its `exp()` returned a
double other than the correctly rounded one for 0.48% of 700,000 arguments. A run on Windows
and a run on Linux can therefore differ in their last bits even within Python.

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

The two CSVs compare exactly, up to CSV number formatting, wherever the two languages share a
maths library. The repository's
[`python/examples/parity_check.py`](https://github.com/pablobernabeu/pilotr/blob/main/python/examples/parity_check.py)
runs this comparison across the worked example designs and reports the maximum difference:
zero for every design with rounded responses, and below 1e-14 for the unrounded continuous
design, where R's CSV writer prints 15 significant digits. At that precision the residual
cannot be told apart from the writer's own rounding, and the script's small tolerance (1e-6)
is there to absorb it. The stricter harness under `tools/parity/` dumps 17 digits and resolves
what is really there. On Windows, it finds a handful of cells one unit in the last place apart,
from `exp()`, and on Linux, where it runs in continuous integration, every cell agrees. That
harness pins the Gaussian cases to recorded hashes, on the observed agreement of AS 241's
`log()` across platforms. The cases are small enough that none has met one of the rare normal
draws that differ. The harness deliberately leaves unpinned the cases whose values pass through
`exp()`, `log()` or `pow()` on the way to the response. Those last bits belong to the maths
library, and so are not the specification's to fix.

`solve_curve` sits a layer above the generator, taking a curve that has already been simulated
and fitting a model to it, and it is held to a different standard on purpose. Its fit calls
`exp` and the normal distribution function at every iteration, so the two engines cannot be
asked to land on the same bits, and a recorded hash would pin whichever maths library happened
to record it. What they are asked for is agreement far closer than any reading of the result
could depend on. The script `tools/parity/solve_cross.py` puts the same fixed curves through
both engines and compares the solved value, its interval, the fitted coefficients and the text
of every refusal. The largest relative difference measured over those cases is 5.0e-15, against
an allowance of 1e-9. On a solved sample size of 220 subjects that allowance comes to two
ten-millionths of a subject. Every refusal message is identical character for character.

The mixed models of the two `power_mixed` functions come from `model_formula` and `model_data`,
which this package ports from the R package. The script `tools/parity/formula_cross.py` gives the
ports the data R simulated for every shipped example and parity case. Both must write the same
formula character for character and the same contrast and product columns bit for bit. The fits
that follow are another matter. `lme4` and `statsmodels` are different programs, the R package
estimates the random-effect correlations that this one leaves out, and the two test with
different statistics. The [power guide](power.md#mixed-effects-power) sets out where their results
diverge.

`power_report` writes a sample-size justification from the counts of a result, the significant
replicates and the replicates behind each rate. It computes each Monte Carlo error with the
functions that compute them in both packages, which return the same doubles for the same counts,
and formats them in the same way. The script `tools/parity/report_cross.py` renders a set of
stored records in both engines and requires the same text, byte for byte.

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
