# Cross-language reproducibility

The point of the shared specification is that one design produces the same data in R and in
Python. The two implementations never call each other but share a numerical contract.

## The shared generator

The R and Python ecosystems use different default random-number generators (R's
Mersenne-Twister with inversion, NumPy's PCG64), so a direct port produces different data from
the same seed. In place of the default generators, `pilotr` provides a shared generator
implemented identically in both languages.

Uniform deviates come from L'Ecuyer's (1988) combined linear congruential generator. Every
intermediate product stays below `2**53`, so the arithmetic is exact in R's IEEE-754 doubles
and in Python integers alike. That holds for seeds within ±(2^53 − 1), the range `validate_spec`
accepts. R reads a JSON number as a double, so beyond that range the two readers would hold
different seeds. Normal deviates use Wichura's (1988) Algorithm AS 241 inverse-CDF, the
algorithm R's `qnorm` uses. The remaining draws (Cholesky-correlated random effects,
inverse-CDF Poisson and ordinal draws, Marsaglia–Tsang gamma draws for the Beta family) are
derived from those two through a documented, identical draw order, which the
[Specification](specification.md) page sets out.

## Where the data can differ

Given the same specification and seed, the two implementations produce identical data wherever
they share a maths library, as on Linux. Continuous integration runs every case of the parity
harness there. The harness holds every case to the last bit except the unrounded dumps of
`crossed_mixed_rt`, `reading_time_continuous` and `stress_exgaussian`, which
`tools/parity/tolerance.json` allows to differ by up to 8 units in the last place, on Linux as on
any other platform.

IEEE-754 fixes the rounding of arithmetic but leaves that of `exp()`, `log()` and `pow()` to the
library. On Windows, Python and R do not share one. CPython takes them from Microsoft's Universal
C Runtime, and R from the MinGW-w64 toolchain it is built with. A value that passes through those
functions can therefore differ in its last bits on Windows, and every family has such a path. The
figures below were measured on Windows 11 with Python 3.13 and R 4.6.1.

- Normal draws pass through `log()` in AS 241's tail branch, which about 15% of draws take.
  About 1.4e-5 of normal draws differed, by 1 to 3 units in the last place (ulp), so Gaussian
  data can differ too.
- Responses of the lognormal families and the ex-Gaussian differed by 1 ulp.
- An unrounded Beta response differed by up to 35 ulp.
- A discrete outcome can in principle move by a whole category, although none did in 50,000
  rows of each discrete family.

One Beta case can differ on Linux too, under R 4.6.0 or later. From that version, R's `^` computes
a cube or a fourth power by repeated multiplication, where Python calls `pow()`, as R did before.
A Beta response whose Gamma shape is exactly 1/3 or 1/4 can therefore differ in its last bits.
Such a shape needs a cell whose rows share one linear predictor, as in a Beta design with no
random effects and a linear predictor of 0, where a `phi` of 0.5 gives shapes of 1/4.

A square is computed by multiplication in both languages. R's `^` computes it that way in every
version, and Python computes a boost whose exponent is exactly 2 as `u * u`. A design with a
`phi` of 1, which gives shapes of 1/2, therefore draws the same data in both. Up to version 0.3.0,
the Python package computed that square with `pow()`, and such a design could differ in its last
bits.

Rounding has removed every continuous difference observed so far, but it guarantees nothing. A
last-bit difference can straddle a rounding boundary, and `round` is not allowed for the discrete
families. The only way to obtain agreement to the last bit in every family is to run both
languages on a platform where they share a maths library, outside the Beta case under R 4.6.0 or
later. The section 'Scope of the guarantee' on the [Specification](specification.md) page gives
the measurements in full.

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
`exp()`, `log()` and `pow()` from the platform's maths library. On Windows, its `exp()` returned a
double other than the correctly rounded one for 0.48% of 700,000 arguments. A run on Windows
and a run on Linux can therefore differ in their last bits even within Python.

Changing the seed changes the draws, while the structure and the ground-truth parameters stay
fixed.

```python exec="true" source="material-block" session="xl"
spec2 = dict(spec, seed=spec["seed"] + 1)
print(simulate(spec).column("score") == simulate(spec2).column("score"))
```

## Verifying parity

Simulate the same specification and seed in both languages and compare. In Python:

```python
from pilotr import simulate, load_spec
simulate(load_spec("design.json")).to_csv("py.csv")
```

In R:

```r
library(pilotr)
write.csv(simulate_design(load_spec("design.json")), "r.csv", row.names = FALSE)
```

The two CSV files agree exactly, up to CSV number formatting, wherever the two languages share a
maths library. The repository's
[`python/examples/parity_check.py`](https://github.com/pablobernabeu/pilotr/blob/main/python/examples/parity_check.py)
compares the CSV files that `python/examples/run_demo.py` and `r/pilotr/examples/run_demo.R`
write to `build/` for six of the worked example designs, all but the ordinal and Poisson ones, so
run both demos first. It reports the maximum difference, which is zero for every design whose
columns are all rounded or whole, and below 1e-14 for `reading_time_continuous`. That design's
continuous predictors are unrounded, and R's CSV writer prints them to 15 significant digits. At
that precision, the residual cannot be told apart from the writer's own rounding, and the
script's tolerance of 1e-6 is there to absorb it.

The stricter harness under `tools/parity/` dumps 17 digits, which shows whether the values
themselves differ. On Windows, it finds a handful of cells one unit in the last place apart, from
`exp()`. On Linux, where it runs in continuous integration, every cell agrees. The harness checks
the Gaussian cases against recorded hashes, relying on the observed agreement of AS 241's `log()`
across platforms. The cases are small enough that none has met one of the rare normal draws that
differ. The cases whose values pass through `exp()`, `log()` or `pow()` on the way to the response
have no recorded hashes, because the maths library determines those last bits.

In the worked crossed reaction-time example, every subject meets every item in both conditions,
so the data have one row per subject, item and condition. The first rows, which the R package
produces too, look like this:

```python exec="true" session="xl"
import sys; sys.path.insert(0, "docs")
from _exec import table, show
```

```python exec="true" source="material-block" session="xl"
from pilotr import simulate, load_spec, pilotr_example

data = simulate(load_spec(pilotr_example("crossed_mixed_rt")))
print(table(data.head(6)))
```

`solve_curve` fits a model to a curve that has already been simulated, and it is held to a
different standard from the generator. Its fit calls `exp` and the normal distribution function
at every iteration. A recorded hash would therefore be valid only for the maths library that
produced it, and the two engines cannot be expected to agree to the last bit. They must instead
agree far more closely than any reading of the result could depend on. The script
`tools/parity/solve_cross.py` puts the same fixed curves through both engines and compares the
solved value, its interval, the fitted coefficients and the text of a selection of error
messages. Each of those error messages is identical character for character. The largest
relative difference measured over those curves is about 5e-15, against an allowance of 1e-9. On a
solved sample size of 220 subjects, that allowance comes to two ten-millionths of a subject.

The mixed models of the two `power_mixed` functions come from `model_formula` and `model_data`,
which this package ports from the R package. The script `tools/parity/formula_cross.py` gives the
ports the data R simulated for every worked example and parity case. Both packages must write the
same formula character for character and the same contrast and product columns bit for bit. The
fits that follow can differ. `lme4` and `statsmodels` are different programs, the R package
estimates the random-effect correlations that this one leaves out, and the two test with different
statistics. The [power guide](power.md#mixed-effects-power) sets out where their results diverge.

`power_report` writes a sample-size justification from a result's counts of significant replicates
and of the replicates behind each rate. It computes each Monte Carlo error with the functions both
packages use to compute these errors, which return the same doubles for the same counts, and
formats the errors in the same way. The script `tools/parity/report_cross.py` renders a set of
stored records in both engines and requires the same text, byte for byte.
