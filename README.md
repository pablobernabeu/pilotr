# simdgp — a cross-language toolkit for simulating experimental & behavioral data

> **Working name `simdgp` (SIMulate Data-Generating Processes) is provisional.**
> Verify availability on CRAN and PyPI before release.

This is the next-generation successor to the *Experimental-data-simulation* Shiny app
(Bernabeu & Lynott, 2020). Where the prototype draws **marginal** distributions only
(`rnorm`/`rbinom`, no effects, no correlations, no random structure), `simdgp` is a
**generative** toolkit built around one idea:

> **One portable design specification → three interchangeable front-ends
> (no-code web app · R package · Python package) that produce *identical* data,
> closing the loop from design → simulate → analyze → simulation-based power
> & design analysis (Type S/M).**

## Why this is novel (the gap, verified against the literature)

| Capability | faux | simstudy | simr / Superpower | **simdgp** |
|---|:--:|:--:|:--:|:--:|
| Generative IV→DV effect sizes | ✗ | ✓ | (from fitted model) | ✓ |
| Crossed by-subject **and** by-item random slopes | ✓ | partial | ✓ | ✓ |
| Realistic distributions (RT/count/ordinal) | ✓ | ✓ | family-dependent | ✓ |
| Simulation-based power + **Type S/M** | ✗ | ✗ | ✓ (power only) | ✓ |
| No-code GUI | ✗ | ✗ | ANOVA only | ✓ (planned) |
| **Python implementation** | ✗ | ✗ | ✗ | **✓** |
| **R = Python bit-identical from one spec** | — | — | — | **✓** |

No single existing tool spans these. The Python column is empty across the board:
SDV *learns from real data*, pyDOE3 only builds design matrices, Faker makes
placeholders, and statsmodels' power module is analytic classical-test only.

## Repository layout

```
toolkit/
  spec/          The portable design-specification format (the conceptual core)
    SPEC.md          human-readable format documentation + the RNG contract
    design.schema.json   JSON Schema for validation
    examples/        worked design specs (between-groups; crossed mixed-effects RT)
  python/        simdgp Python package (runnable; pure-Python generative core)
  r/             simdgp R package (mirrors the Python core exactly)
  docs/
    positioning.md   the Behavior Research Methods positioning statement / abstract
```

## Cross-language reproducibility (the hard part, solved)

Native RNGs differ across ecosystems (R uses Mersenne-Twister + inversion; NumPy uses
PCG64), so naive ports never match. `simdgp` instead ships a **shared generator**
implemented identically in both languages:

* **Uniforms:** L'Ecuyer (1988) combined linear congruential generator. All integer
  arithmetic stays below 2^53, so it is *exact* in IEEE-754 doubles (R) and Python ints.
* **Normals:** Wichura's (1988) Algorithm AS 241 inverse-CDF — the same algorithm R's
  `qnorm()` uses — so deviates agree to full double precision.
* **Everything else** (Cholesky for correlated random effects, Poisson/Bernoulli/ordinal)
  is derived from those two via inverse-CDF transforms, in a **documented draw order**
  (see `spec/SPEC.md`). Same spec + same seed ⇒ identical data in R and Python.

This is the v0.x research-grade engine (auditable, dependency-free). The production
roadmap upgrades the generator to a counter-based RNG (Philox/Threefry) via a compiled
backend for speed and parallel streams, preserving the same cross-language contract.

## Quick start

```bash
# Python: simulate both designs + classical simulation-based power (Type S/M)
python toolkit/python/examples/run_demo.py

# R: the same, bit-for-bit
Rscript toolkit/r/simdgp/examples/run_demo.R

# Prove R and Python produce identical data (max abs diff = 0):
python toolkit/python/examples/parity_check.py

# R: crossed mixed-effects simulation-based power via lme4/lmerTest
Rscript toolkit/r/simdgp/examples/run_power_mixed.R

# R: validate the generative model -- maximal lmer recovers the specified parameters
Rscript toolkit/r/simdgp/examples/validate_recovery.R

# Python validation suite
python -m pytest toolkit/python/tests -q

# Realistic distributions: ordinal (Likert) + Poisson counts
python toolkit/python/examples/families_demo.py

# Power-vs-N curves + the publication figure
python toolkit/python/examples/power_curves.py        # Gaussian curve  -> build/*.csv
Rscript toolkit/r/simdgp/examples/power_curve_mixed.R # crossed mixed curve (slow, lme4)
Rscript toolkit/r/simdgp/examples/plot_power_curves.R # -> build/power_curves.png

# faux equivalence (the "why not faux?" rebuttal)
Rscript toolkit/r/simdgp/examples/equivalence_faux.R
```

> **R↔Python asymmetry (by design, for now).** Data *generation* is identical across
> languages (proven). For *analysis*, crossed mixed-effects power uses R's mature
> `lme4`/`lmerTest`; the Python analysis backends (`statsmodels` MixedLM is single-grouping
> only — `pymer4`/`formulaic` are the route to crossed models) are on the roadmap. The
> two-group Gaussian power backend runs identically in both languages.
