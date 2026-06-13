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
    inst/app/      the no-code Shiny app, bundled in the package (simdgp::run_app())
  app-lite/      serverless (shinylive/webR) build of the light path -> static site
  docs/
    positioning.md   the Behavior Research Methods positioning statement / abstract
```

## One model, three interfaces

The same design spec drives a **no-code web app**, an **R package**, and a **Python
package**. The app is a thin client: every control writes into the portable JSON spec,
which you can download and run unchanged in either package to get identical data.

```r
# launch the no-code app locally (installed package)
simdgp::run_app()
# ...or from source:
shiny::runApp("toolkit/r/simdgp/inst/app")
```

## Deployment & concurrency

R is single-threaded: one R process runs one computation at a time, and a heavy
simulation-based **power** run (hundreds–thousands of model fits) blocks every other user
sharing that process. So the architecture is deliberately split by how the tool is used:

| Path | How | Concurrency | Use for |
|---|---|---|---|
| **Installable** (primary) | `simdgp::run_app()` / `python -m simdgp`; CRAN/PyPI | **unbounded** — each user, own machine & cores | real work, esp. heavy power (parallelise across cores) |
| **Serverless demo** | `app-lite/` exported with shinylive → static site on GitHub Pages | **unbounded** — each browser computes (WebAssembly) | zero-cost "try it now" link: design + simulate + Gaussian power |
| Shared hosted instance | shinyapps.io / ShinyProxy | low / costly | avoid as the main channel (blocking; the prototype's 25 hr-month free tier) |

The installable app runs power **asynchronously** (via `future`/`promises`) so it never
blocks; the serverless build is single-user-per-browser, so it runs synchronously. Both are
driven by the *same* spec, so a user can design in the browser demo and run heavy power
locally from the downloaded spec.

```bash
# build the serverless static site (downloads webR assets on first run)
Rscript toolkit/app-lite/build_shinylive.R   # -> toolkit/build/shinylive-demo/ (deploy to gh-pages)
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

# competitor equivalence (the "why not faux / simstudy?" rebuttals)
Rscript toolkit/r/simdgp/examples/equivalence_faux.R
Rscript toolkit/r/simdgp/examples/equivalence_simstudy.R

# crossed mixed-effects power in Python (statsmodels) -- R/Python capability parity
python toolkit/python/examples/power_mixed_demo.py
```

> **R↔Python coverage.** Data *generation* is bit-identical across languages (proven). Both
> ecosystems also run crossed mixed-effects power from the same spec: R via `lme4`/`lmerTest`
> (REML, correlated random effects) and Python via `statsmodels` MixedLM (ML, crossed
> variance components). The two backends are *close but not identical* by construction
> (different estimator/optimizer); the two-group Gaussian power backend is identical in both.
