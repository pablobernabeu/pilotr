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
[`qnorm()`](https://rdrr.io/r/stats/Normal.html) uses, so deviates agree
to full double precision. The remaining draws (Cholesky-correlated
random effects, inverse-CDF Poisson and ordinal draws, Marsaglia–Tsang
gamma draws for the Beta family, a Fisher–Yates item sampler for partial
crossing) are derived from those two through a documented, identical
consumption order.

As a result, given the same specification and seed, the R and Python
implementations produce identical data sets: bit for bit for rounded and
Gaussian responses, and within a few units in the last place where an
unrounded family applies [`exp()`](https://rdrr.io/r/base/Log.html) or
[`log()`](https://rdrr.io/r/base/Log.html), as the parity check below
quantifies. The specification is the source of truth, and every
implementation answers to it.

## Determinism in R

The same specification and seed always give the same data.

``` r

spec <- build_spec(list(
  name = "demo", seed = 2024, design_kind = "between", n_subject = 200,
  factor_name = "group", lev1 = "a", lev2 = "b",
  intercept = 100, effect = 5, family = "gaussian",
  resp_name = "score", sigma = 10))

isTRUE(all.equal(simulate_design(spec), simulate_design(spec)))
```

    [1] TRUE

Changing the seed changes the draws, while the structure and the
ground-truth parameters stay fixed.

``` r

spec2 <- spec; spec2$seed <- spec$seed + 1L
identical(simulate_design(spec)$score, simulate_design(spec2)$score)
```

    [1] FALSE

## The same design in Python

The Python package reads the identical specification and produces the
identical data. The code below is not evaluated here. It shows the
Python you would run, and its output matches the R data to the last bit.

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
check’s 1e-6 tolerance is there to absorb. The stricter harness dumps 17
digits and resolves what is really there: a handful of cells apart by a
few units in the last place, from
[`exp()`](https://rdrr.io/r/base/Log.html), whose rounding IEEE-754 does
not fix. There are no categorical mismatches anywhere. The uniform
stream is exact by construction, since every intermediate product stays
below $`2^{53}`$, and the harness pins every path built from IEEE-exact
arithmetic alone to a recorded hash in `tools/parity/golden.json`. The
cases whose values pass through
[`exp()`](https://rdrr.io/r/base/Log.html) or
[`log()`](https://rdrr.io/r/base/Log.html) are deliberately left
unpinned, because those last bits belong to the maths library and so are
not the specification’s to fix. The cross-language comparison still
gates them on every run. The parity check can be rerun on any new
platform, so a design verified on a laptop can be confirmed on an HPC
node before a large run.

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

## Why this matters

Because the specification is the source of truth and the generator is
shared, a design can be authored once (in the app, in R or in Python)
and reproduced everywhere. It can be preregistered as a small JSON file
that anyone can run, scaled from a laptop to a cluster without drift and
handed between collaborators who work in different languages.

This property turns an informal request to simulate some data into a
reproducible, shareable methodological artefact.
