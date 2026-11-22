# Cross-language reproducibility

``` r

library(pilotr)
```

## The challenge

The R and Python ecosystems use different default random-number
generators (R’s Mersenne-Twister with inversion, NumPy’s PCG64), so a
direct port produces different data from the same seed. That would
undermine the guarantee that a design authored once can be reproduced
anywhere.

## A shared generator

In place of the default generators, pilotr provides a shared generator
implemented identically in both languages. Uniform deviates come from
L’Ecuyer’s (1988) combined linear congruential generator. Every
intermediate product stays below $`2^{53}`$, so the arithmetic is exact
in R’s IEEE-754 doubles and in Python integers alike. Normal deviates
use Wichura’s (1988) Algorithm AS 241 inverse-CDF, the algorithm R’s
[`qnorm()`](https://rdrr.io/r/stats/Normal.html) uses. The remaining
draws (Cholesky-correlated random effects, inverse-CDF Poisson and
ordinal draws, Marsaglia–Tsang gamma draws for the Beta family, a
Fisher–Yates item sampler for partial crossing) are derived from those
two through a documented, identical draw order.

As a result, given the same specification and seed, the R and Python
implementations produce identical data wherever they share a maths
library, as on Linux, except in a rare Beta case under R 4.6.0 or later.
Continuous integration runs every case of the parity harness on Linux.
The harness holds every case to the last bit except the unrounded dumps
of `crossed_mixed_rt`, `reading_time_continuous` and
`stress_exgaussian`, which `tools/parity/tolerance.json` allows to
differ by up to 8 units in the last place, on Linux as on any other
platform.

IEEE-754 fixes the rounding of arithmetic but leaves that of
[`exp()`](https://rdrr.io/r/base/Log.html),
[`log()`](https://rdrr.io/r/base/Log.html) and `pow()` to the library.
On Windows, R and Python do not share one. There, a value that passes
through those functions can differ in its last bits, and every family
has such a path. The figures below were measured on Windows 11 with R
4.6.1 and Python 3.13.

- Normal draws pass through [`log()`](https://rdrr.io/r/base/Log.html)
  in AS 241’s tail branch, which about 15% of draws take. About 1.4e-5
  of normal draws differed, by 1 to 3 units in the last place (ulp), so
  Gaussian data can differ too.
- Responses of the lognormal families and the ex-Gaussian differed by 1
  ulp.
- An unrounded Beta response differed by up to 35 ulp.
- A discrete outcome can in principle move by a whole category, although
  none did in 50,000 rows of each discrete family.

One Beta case can differ on Linux too, under R 4.6.0 or later. From that
version, R’s `^` computes a cube or a fourth power by repeated
multiplication, where Python calls `pow()`, as R did before. A Beta
response whose Gamma shape is exactly 1/3 or 1/4 can therefore differ in
its last bits. Such a shape needs a cell whose rows share one linear
predictor, as in a Beta design with no random effects and a linear
predictor of 0, where a `phi` of 0.5 gives shapes of 1/4.

A square is computed by multiplication in both languages. R’s `^`
computes it that way in every version, and Python computes a boost whose
exponent is exactly 2 as `u * u`. A design with a `phi` of 1, which
gives shapes of 1/2, therefore draws the same data in both. Up to
version 0.3.0, the Python package computed that square with `pow()`, and
such a design could differ in its last bits.

Rounding has removed every continuous difference observed so far, but it
guarantees nothing. A last-bit difference can straddle a rounding
boundary, and `round` is not allowed for the discrete families.
[`build_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/build_spec.md),
which builds the no-code app’s specifications, rounds the Gaussian,
lognormal, shifted lognormal and ex-Gaussian responses to four decimals
and leaves a Beta response unrounded, for the reason
[`?build_spec`](https://pablobernabeu.github.io/pilotr/r/reference/build_spec.md)
gives. The only way to obtain agreement to the last bit in every family
is to run both languages on a platform where they share a maths library,
outside the Beta case under R 4.6.0 or later. The section [Scope of the
guarantee](https://github.com/pablobernabeu/pilotr/blob/main/spec/SPEC.md#scope-of-the-guarantee)
of the specification gives the measurements in full.

## Determinism in R

On one platform, the same specification and seed always give the same
data.

``` r

spec <- build_spec(list(
  name = "demo", seed = 2024, design_kind = "between", n_subject = 200,
  factor_name = "group", lev1 = "a", lev2 = "b",
  intercept = 100, effect = 5, family = "gaussian",
  resp_name = "score", sigma = 10))

identical(simulate_design(spec), simulate_design(spec))
```

    [1] TRUE

Across operating systems, the same holds on the terms set out above,
since R takes [`exp()`](https://rdrr.io/r/base/Log.html),
[`log()`](https://rdrr.io/r/base/Log.html) and `pow()` from the
platform’s maths library, as Python does.

Changing the seed changes the draws, while the structure and the
ground-truth parameters stay fixed.

``` r

spec2 <- spec; spec2$seed <- spec$seed + 1L
identical(simulate_design(spec)$score, simulate_design(spec2)$score)
```

    [1] FALSE

## The same design in Python

The Python package reads the same specification. The code below, which
is not evaluated here, shows the Python you would run. Its output
matches the R data to the last bit wherever the two share a maths
library.

``` python
from pilotr import simulate, load_spec

# the same design.json the R package (and the app) use
data = simulate("design.json")
```

The repository includes a parity check that compares the CSV files the R
and Python demo scripts write for six of the worked example designs, all
but the ordinal and Poisson ones. It reports the maximum absolute
difference per design, reads values at the precision a written file
carries and has a tolerance of 1e-6. The table below summarises it.

| design     | rows | columns | max_abs_diff | categorical_mismatches |
|:-----------|-----:|--------:|:-------------|-----------------------:|
| between    |   64 |       3 | 0            |                      0 |
| crossed    | 1440 |       4 | 0            |                      0 |
| continuous | 4000 |       7 | 4.885e-15    |                      0 |
| nested     | 4800 |       5 | 0            |                      0 |
| beta       |  200 |       3 | 0            |                      0 |
| partial    | 1440 |       4 | 0            |                      0 |

The difference is exactly zero for five of the six designs. The sixth,
`continuous` (`reading_time_continuous`), has a rounded response that
agrees exactly. Its continuous predictors are unrounded, and R’s CSV
writer prints them to 15 significant digits. Those predictors differ by
about 5e-15. At that precision, the residual cannot be told apart from
the writer’s own rounding, and the tolerance of 1e-6 is there to absorb
it. The parity check can be rerun on any new platform, so a design
verified on a laptop can be confirmed on an HPC node before a large run.

The bit-identical claim rests on a stricter check, the harness under
`tools/parity/`. It dumps every value at 17 significant digits, which
shows whether the values themselves differ. It accounts for each
differing cell in units in the last place against a recorded allowance
and must pass in continuous integration. On Windows, it finds a handful
of cells one unit in the last place apart, from
[`exp()`](https://rdrr.io/r/base/Log.html), and no categorical
mismatches. On Linux, where it runs in continuous integration, every
cell agrees. The uniform stream is exact by construction. The harness
checks the Gaussian cases against recorded hashes in
`tools/parity/golden.json`, relying on the observed agreement of AS
241’s [`log()`](https://rdrr.io/r/base/Log.html) across platforms. The
cases are small enough that none has met one of the rare normal draws
that differ. The cases whose values pass through
[`exp()`](https://rdrr.io/r/base/Log.html),
[`log()`](https://rdrr.io/r/base/Log.html) or `pow()` on the way to the
response have no recorded hashes, because the maths library determines
those last bits. They must still pass the cross-language comparison on
every run.

## The analysis layer

[`solve_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/solve_curve.md)
fits a model to a curve that has already been simulated, and it is held
to a different standard from the generator. Its fit calls
[`exp()`](https://rdrr.io/r/base/Log.html) and the normal distribution
function at every iteration. A recorded hash would therefore be valid
only for the maths library that produced it, and the two engines cannot
be expected to agree to the last bit. They must instead agree far more
closely than any reading of the result could depend on. The script
`tools/parity/solve_cross.py` puts the same fixed curves through both
engines and compares the solved value, its interval, the fitted
coefficients and the text of a selection of error messages. Each of
those error messages is identical character for character. The largest
relative difference measured over those curves is about 5e-15, against
an allowance of 1e-9. On a solved sample size of 220 subjects, that
allowance comes to two ten-millionths of a subject.

The mixed models of the two
[`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md)
functions come from
[`model_formula()`](https://pablobernabeu.github.io/pilotr/r/reference/model_formula.md)
and
[`model_data()`](https://pablobernabeu.github.io/pilotr/r/reference/model_data.md),
which the Python package ports from this one. The script
`tools/parity/formula_cross.py` gives the Python ports the data R
simulated for every worked example and parity case. Both packages must
write the same formula character for character and the same contrast and
product columns bit for bit. The fits that follow can differ. `lme4` and
`statsmodels` are different programs, this package estimates the
random-effect correlations that the Python one leaves out, and the two
test with different statistics. The Python power guide sets out where
their results diverge.

[`power_report()`](https://pablobernabeu.github.io/pilotr/r/reference/power_report.md)
writes a sample-size justification from a result’s counts of significant
replicates and of the replicates behind each rate. It computes each
Monte Carlo error with the functions both packages use to compute these
errors, which return the same doubles for the same counts, and formats
the errors in the same way. The script `tools/parity/report_cross.py`
renders a set of stored records in both engines and requires the same
text, byte for byte.

## Why this matters

Because every implementation must follow the specification and the
generator is shared, a design authored once in the app, in R or in
Python can be reproduced everywhere. It can be preregistered as a small
JSON file that anyone can run, scaled from a laptop to a cluster and
handed between collaborators who work in different languages. Its data
agree to the last bit wherever the machines share a maths library,
outside the Beta case under R 4.6.0 or later.
