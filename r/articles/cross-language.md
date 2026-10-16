# Cross-language reproducibility

``` r

library(pilotr)
```

## The challenge

A toolkit that spans R and Python faces a subtle obstacle. The two
ecosystems use different default random-number generators (R’s
Mersenne-Twister with inversion, NumPy’s PCG64), so a direct port
produces different data from the same seed. That would undermine the
guarantee that a design authored once can be reproduced anywhere.

## A shared generator

pilotr sets aside the native RNGs and provides a shared generator
implemented identically in both languages. Uniform deviates come from
L’Ecuyer’s (1988) combined linear congruential generator, whose
arithmetic stays below $`2^{53}`$ and is therefore exact in IEEE-754
doubles (R) and Python integers alike. Normal deviates use Wichura’s
(1988) Algorithm AS 241 inverse-CDF, the algorithm R’s
[`qnorm()`](https://rdrr.io/r/stats/Normal.html) uses. The remaining
draws (Cholesky-correlated random effects, inverse-CDF Poisson and
ordinal draws, Marsaglia–Tsang gamma draws for the Beta family, a
Fisher–Yates item sampler for partial crossing) are derived from those
two through a documented, identical consumption order.

As a result, given the same specification and seed, the R and Python
implementations produce identical data sets wherever they share a maths
library, as on Linux. Continuous integration checks every case of the
parity harness to the last bit there. IEEE-754 fixes the rounding of
arithmetic but leaves that of
[`exp()`](https://rdrr.io/r/base/Log.html),
[`log()`](https://rdrr.io/r/base/Log.html) and `pow()` to the library,
and on Windows R and Python do not share one. There, a value that passes
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

One Beta case can differ on Linux too. From R 4.6.0, R’s `^` computes a
cube or a fourth power by repeated multiplication, where Python calls
`pow()`. A Beta response whose Gamma shape is exactly 1/3 or 1/4 can
therefore differ in its last bits on any platform. Such a shape needs a
cell whose rows share one linear predictor, as in a Beta design with no
random effects, a linear predictor of 0 and `phi` set to 0.5.

Rounding has removed every continuous difference observed so far, but it
guarantees nothing, since a last-bit difference can straddle a rounding
boundary, and `round` is refused for the discrete families.
[`build_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/build_spec.md),
which builds the no-code app’s specifications, rounds the Gaussian and
lognormal responses to four decimals and leaves a Beta response
unrounded, for the reason
[`?build_spec`](https://pablobernabeu.github.io/pilotr/r/reference/build_spec.md)
gives. Agreement to the last bit in every family has one route: run both
languages on a platform where they share a maths library, outside that
Beta case. The specification is the source of truth, and every
implementation answers to it. Its section [Scope of the
guarantee](https://github.com/pablobernabeu/pilotr/blob/main/spec/SPEC.md#scope-of-the-guarantee)
gives the measurements in full.

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

The Python package reads the same specification and produces the same
data. The code below is not evaluated here. It shows the Python you
would run, and its output matches the R data to the last bit wherever
the two share a maths library.

``` python
from pilotr import simulate, load_spec

# the same design.json the R package (and the app) use
data = simulate("design.json")
```

The repository includes a parity check that simulates the worked example
designs in both languages and reports the maximum absolute difference
per design. This is the CSV-level check, read at the precision a written
file carries and with a tolerance of 1e-6. The stricter one is the
harness under `tools/parity/`, which dumps every value at 17 significant
digits, accounts for each differing cell in units in the last place
against a recorded allowance, and runs as a gate in continuous
integration. The table below is the readable summary, and that harness
is what the bit-identical claim rests on.

| design     | rows | columns | max_abs_diff | categorical_mismatches |
|:-----------|-----:|--------:|:-------------|-----------------------:|
| between    |   64 |       3 | 0            |                      0 |
| crossed    | 1440 |       4 | 0            |                      0 |
| continuous | 4000 |       7 | 4.885e-15    |                      0 |
| nested     | 4800 |       5 | 0            |                      0 |
| beta       |  200 |       3 | 0            |                      0 |
| partial    | 1440 |       4 | 0            |                      0 |

The difference is exactly zero for every design whose responses are
rounded, and about 5e-15 for the unrounded continuous design, where R’s
CSV writer prints 15 significant digits. At that precision the residual
cannot be told apart from the writer’s own rounding, and it is what the
1e-6 tolerance of the check is there to absorb. The stricter harness
dumps 17 digits and resolves what is really there. On Windows, it finds
a handful of cells one unit in the last place apart, from
[`exp()`](https://rdrr.io/r/base/Log.html), and no categorical
mismatches. On Linux, where it runs in continuous integration, every
cell agrees. The uniform stream is exact by construction, since every
intermediate product stays below $`2^{53}`$. The harness pins the
Gaussian cases to recorded hashes in `tools/parity/golden.json`, on the
observed agreement of AS 241’s
[`log()`](https://rdrr.io/r/base/Log.html) across platforms. The cases
are small enough that none has met one of the rare normal draws that
differ. The cases whose values pass through
[`exp()`](https://rdrr.io/r/base/Log.html),
[`log()`](https://rdrr.io/r/base/Log.html) or `pow()` on the way to the
response are deliberately left unpinned, because those last bits belong
to the maths library and so are not the specification’s to fix. The
cross-language comparison still gates them on every run. The parity
check can be rerun on any new platform, so a design verified on a laptop
can be confirmed on an HPC node before a large run.

## The analysis layer

Everything above concerns the data the generator produces.
[`solve_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/solve_curve.md)
sits a layer above it, taking a curve that has already been simulated
and fitting a model to it, and it is held to a different standard on
purpose. Its fit calls [`exp()`](https://rdrr.io/r/base/Log.html) and
the normal distribution function at every iteration, so the two engines
cannot be asked to land on the same bits, and a recorded hash would pin
whichever maths library happened to record it. What they are asked for
is agreement far closer than any reading of the result could depend on.
The script `tools/parity/solve_cross.py` puts the same fixed curves
through both engines and compares the solved value, its interval, the
fitted coefficients and the text of every refusal. The largest relative
difference measured over those cases is 5.0e-15, against an allowance of
1e-9. On a solved sample size of 220 subjects that allowance comes to
two ten-millionths of a subject. Every refusal message is identical
character for character.

The mixed models of the two
[`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md)
functions come from
[`model_formula()`](https://pablobernabeu.github.io/pilotr/r/reference/model_formula.md)
and
[`model_data()`](https://pablobernabeu.github.io/pilotr/r/reference/model_data.md),
which the Python package ports from this one. The script
`tools/parity/formula_cross.py` gives the Python ports the data R
simulated for every shipped example and parity case. Both must write the
same formula character for character and the same contrast and product
columns bit for bit. The fits that follow are another matter. `lme4` and
`statsmodels` are different programs, this package estimates the
random-effect correlations that the Python one leaves out, and the two
test with different statistics. The Python power guide sets out where
their results diverge.

[`power_report()`](https://pablobernabeu.github.io/pilotr/r/reference/power_report.md)
writes a sample-size justification from the counts of a result, the
significant replicates and the replicates behind each rate. It computes
each Monte Carlo error with the functions that compute them in both
packages, which return the same doubles for the same counts, and formats
them in the same way. The script `tools/parity/report_cross.py` renders
a set of stored records in both engines and requires the same text, byte
for byte.

## Why this matters

Because the specification is the source of truth and the generator is
shared, a design can be authored once (in the app, in R or in Python)
and reproduced everywhere. It can be preregistered as a small JSON file
that anyone can run, scaled from a laptop to a cluster and handed
between collaborators who work in different languages, with data that
agree to the last bit wherever the machines share a maths library.

This property turns an informal request to simulate some data into a
reproducible, shareable methodological artefact.
