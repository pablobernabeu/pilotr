# pilotr, a cross-language toolkit for simulating experimental and behavioral data

`pilotr` lets researchers pilot a study before they run it.

The package is a successor to the *Experimental-data-simulation* Shiny app
(Bernabeu & Lynott, 2020), archived on Zenodo (doi:10.5281/zenodo.10615953). The
prototype draws marginal distributions only (`rnorm`/`rbinom`, with no effects, no
correlations, and no random structure). `pilotr` is a generative toolkit organised
around a single idea. One portable design specification drives three interchangeable
front-ends, a no-code web application, an R package, and a Python package, that produce
identical data. Together they close the loop from design to simulation to analysis, and
on to simulation-based power and design analysis (Type S/M).

## Position relative to existing tools (verified against the literature)

| Capability | faux | simstudy | simr / Superpower | **pilotr** |
|---|:--:|:--:|:--:|:--:|
| Generative IV→DV effect sizes | ✗ | ✓ | (from fitted model) | ✓ |
| Crossed by-subject **and** by-item random slopes | ✓ | partial | ✓ | ✓ |
| Realistic distributions (RT/count/ordinal) | ✓ | ✓ | family-dependent | ✓ |
| Simulation-based power + **Type S/M** | ✗ | ✗ | ✓ (power only) | ✓ |
| No-code GUI | ✗ | ✗ | ANOVA only | ✓ |
| **Python implementation** | ✗ | ✗ | ✗ | **✓** |
| **R = Python bit-identical from one spec** | — | — | — | **✓** |

No single existing tool spans these capabilities. The Python column is empty across the
board. SDV learns from real data, pyDOE3 builds design matrices, Faker produces
placeholder values, and the statsmodels power module covers analytic classical tests only.

## Repository layout

```
pilotr/  (repo root)
  spec/          The portable design-specification format (the conceptual core)
    SPEC.md          human-readable format documentation + the RNG contract
    design.schema.json   JSON Schema for validation
    examples/        worked design specs (between-groups; crossed mixed-effects RT)
  python/        pilotr Python package (runnable; pure-Python generative core)
  r/             pilotr R package (mirrors the Python core exactly)
    inst/app/      the no-code Shiny app, bundled in the package (pilotr::run_app())
  app-lite/      serverless (shinylive/webR) build of the light path -> static site
  docs/
    positioning.md   the Behavior Research Methods positioning statement / abstract
```

## One model, three interfaces

The same design spec drives a no-code web app, an R package, and a Python package. The
app is a thin client. Every control writes into the portable JSON spec, which can be
downloaded and run unchanged in either package to obtain identical data.

```r
# launch the no-code app locally (installed package)
pilotr::run_app()
# ...or from source:
shiny::runApp("r/pilotr/inst/app")
```

## Deployment and concurrency

R is single-threaded. One R process runs one computation at a time, and a heavy
simulation-based power run (hundreds to thousands of model fits) blocks every other user
sharing that process. The architecture is therefore split according to how the tool is
used.

| Path | How | Concurrency | Use for |
|---|---|---|---|
| Installable (primary) | `pilotr::run_app()` in R, or `import pilotr` for Python scripting (installed from source until released on CRAN and PyPI) | unbounded, each user on their own machine and cores | real work, especially heavy power runs parallelised across cores |
| Serverless demo | `app-lite/` exported with shinylive to a static site on GitHub Pages | unbounded, each browser computes via WebAssembly | a low-cost link for design, simulation, and Gaussian power |
| Shared hosted instance | shinyapps.io or ShinyProxy | low and costly | best avoided as the main channel, since it blocks and the prototype's free tier allowed 25 hours per month |

The installable app runs power asynchronously (via `future`/`promises`) so that it does
not block. The serverless build is single-user-per-browser and runs synchronously. Both
are driven by the same spec, so a user can design in the browser demo and run heavy power
locally from the downloaded spec.

```bash
# build the serverless static site (downloads webR assets on first run)
Rscript app-lite/build_shinylive.R   # -> build/shinylive-demo/

# The published website is built by CI (.github/workflows/site.yml): the pkgdown docs at
# https://pablobernabeu.github.io/pilotr/ and this demo at https://pablobernabeu.github.io/pilotr/demo/
```

## Running at scale (HPC / SLURM)

Simulation-based power and precision analyses are embarrassingly parallel, so they scale
well on a cluster. `hpc/` holds a robust SLURM array job (`precision_array.slurm` and its
runner `precision_array.R`). It runs one task per sample size, with replicates
parallelised across cores via `mclapply` and results written to project storage. A
reference deployment on the Oxford ARC cluster is organised as follows.

- **Project (code/scripts):** `~/pilotr_toolkit/` in home. **Heavy material (R library +
  results):** `/data/<project>/pilotr_toolkit/{Rlib,results,logs}` in project storage, since
  home quota is small.
- A one-time bootstrap installs `lme4` and `lmerTest` into the data-area library
  (`R_LIBS`). The R module already provides `jsonlite`.
- Submit the sweep: `sbatch hpc/precision_array.slurm`
  (smoke test: `sbatch --export=ALL,N_SIMS=4 --array=0 --partition=devel precision_array.slurm`).
- Each task writes one `precision_N<n>.csv`. These combine into a full precision-vs-N curve
  at a resolution far beyond a laptop.

The simulation core is bit-identical across machines and R versions. ARC reproduces local
output exactly, so HPC runs are reproducible from the same spec and seed.

## Cross-language reproducibility

Native RNGs differ across ecosystems. R uses Mersenne-Twister with inversion, whereas
NumPy uses PCG64, so naive ports never match. `pilotr` instead ships a shared generator
implemented identically in both languages.

* **Uniforms:** L'Ecuyer (1988) combined linear congruential generator. All integer
  arithmetic stays below 2^53, so it is exact in IEEE-754 doubles (R) and Python ints.
* **Normals:** Wichura's (1988) Algorithm AS 241 inverse-CDF, the same algorithm R's
  `qnorm()` uses, so deviates agree to full double precision.
* **Everything else** (Cholesky for correlated random effects, Poisson/Bernoulli/ordinal)
  is derived from those two via inverse-CDF transforms, in a documented draw order
  (see `spec/SPEC.md`). The same spec and seed yield identical data in R and Python.

This is the v0.x research-grade engine, which is auditable and dependency-free. The
production roadmap upgrades the generator to a counter-based RNG (Philox/Threefry) via a
compiled backend for speed and parallel streams, preserving the same cross-language
contract.

## Quick start

```bash
# Python: simulate both designs + classical simulation-based power (Type S/M)
python python/examples/run_demo.py

# R: the same, bit-for-bit
Rscript r/pilotr/examples/run_demo.R

# Check that R and Python produce identical data (max abs diff = 0)
python python/examples/parity_check.py

# R: crossed mixed-effects simulation-based power via lme4/lmerTest
Rscript r/pilotr/examples/run_power_mixed.R

# R: validate the generative model (a maximal lmer fit recovers the specified parameters)
Rscript r/pilotr/examples/validate_recovery.R

# Python validation suite
python -m pytest python/tests -q

# Realistic distributions: ordinal (Likert) + Poisson counts
python python/examples/families_demo.py

# Power-vs-N curves + the publication figure
python python/examples/power_curves.py        # Gaussian curve  -> build/*.csv
Rscript r/pilotr/examples/power_curve_mixed.R # crossed mixed curve (slow, lme4)
Rscript r/pilotr/examples/plot_power_curves.R # -> build/power_curves.png

# equivalence with faux and simstudy where they overlap
Rscript r/pilotr/examples/equivalence_faux.R
Rscript r/pilotr/examples/equivalence_simstudy.R

# crossed mixed-effects power in Python (statsmodels), for R and Python capability parity
python python/examples/power_mixed_demo.py

# continuous predictors, interactions, and continuous random slopes, with a precision-based
# ROPE design analysis, an N-sweep, and a brms bridge (see docs/mixed_models_and_design_analysis.md)
Rscript r/pilotr/examples/precision_design_analysis.R
```

> **R↔Python coverage.** Data generation is bit-identical across languages (proven), and
> both ecosystems run crossed mixed-effects power from the same spec. The LMM estimators,
> however, differ. R (`lme4`/`lmerTest`, REML, correlated random effects) is the reference.
> Python (`statsmodels` MixedLM, crossed variance components) overstates random-slope
> variance, so it is conservative for random-slope designs (for the crossed design, power
> about 0.48 versus about 0.73 from lme4), although it recovers fixed effects correctly
> (mean estimate about 0.048 versus 0.05). The
> two-group Gaussian power backend is identical in both languages. For correlated random
> slopes, R is the recommended choice.
