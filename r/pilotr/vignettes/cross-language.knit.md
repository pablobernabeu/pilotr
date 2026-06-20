---
title: "Cross-language reproducibility"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Cross-language reproducibility}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




``` r
library(pilotr)
```

## The problem

A toolkit that spans R and Python faces a subtle obstacle: the two ecosystems use different
default random-number generators (R's Mersenne-Twister with inversion; NumPy's PCG64), so a
naive port produces *different* data from the same seed. That breaks the promise that a
design authored once can be reproduced anywhere.

## The solution: a shared generator

pilotr sidesteps native RNGs entirely and ships a **shared generator** implemented identically
in both languages:

- **Uniform deviates** come from L'Ecuyer's (1988) combined linear congruential generator,
  whose arithmetic stays below $2^{53}$ and is therefore exact in IEEE-754 doubles (R) and
  Python integers alike.
- **Normal deviates** use Wichura's (1988) Algorithm AS 241 inverse-CDF — the algorithm R's
  `qnorm()` uses — so deviates agree to full double precision.
- **Everything else** (Cholesky-correlated random effects, inverse-CDF Poisson and ordinal
  draws, Marsaglia–Tsang gamma draws for the Beta family, a Fisher–Yates item sampler for
  partial crossing) is derived from those two through a **documented, identical consumption
  order**.

The consequence: given the same specification and seed, the R and Python implementations
produce **bit-identical** data sets. The *specification*, not the implementation, is the
source of truth.

## Determinism in R

The same spec and seed always give the same data:


``` r
spec <- build_spec(list(
  name = "demo", seed = 2024, design_kind = "between", n_subject = 200,
  factor_name = "group", lev1 = "a", lev2 = "b",
  intercept = 100, effect = 5, family = "gaussian", resp_name = "score", sigma = 10))

isTRUE(all.equal(simulate_design(spec), simulate_design(spec)))
#> [1] TRUE
```

Change the seed and the draws change; the structure (and the ground-truth parameters) stay
fixed:


``` r
spec2 <- spec; spec2$seed <- spec$seed + 1L
identical(simulate_design(spec)$score, simulate_design(spec2)$score)
#> [1] FALSE
```

## The same design in Python

The Python package reads the identical specification and produces the identical data. Nothing
below is evaluated here — it is the Python you would run — but the output matches the R data
to the last bit:

```python
from pilotr import simulate, load_spec

# the same design.json the R package (and the app) use
data = simulate("design.json")
```

The repository ships a parity check that simulates both designs in both languages and reports
the maximum absolute difference — **0.0** across thousands of rows, with no categorical
mismatches. The simulation core is also bit-identical *across machines and R versions*, so a
design verified on a laptop reproduces exactly on an HPC node.

## Why this matters

Because the specification is the source of truth and the generator is shared, a design can be:

- **authored once** (in the app, in R, or in Python) and reproduced everywhere;
- **preregistered** as a small JSON file that anyone can run;
- **scaled** from a laptop to a cluster without drift;
- **handed between collaborators** who work in different languages.

This is the property that turns "simulate some data" into a reproducible, shareable
methodological artifact.
