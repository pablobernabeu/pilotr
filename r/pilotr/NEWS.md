# pilotr 0.2.1

Released 2026-07-23.

## Documentation

* A package-level help topic, `?pilotr`, introduces the workflow and groups
  every export by the stage it belongs to.
* The reference examples that cut the replicate count so they run quickly now
  say so, and each `n_sims` argument documents the Monte Carlo error and a
  usable recommendation for real planning.
* The `precision_curve()` example passes a region of practical equivalence
  clearly inside the effect it declares. It previously fell back to a default
  equal to that effect, the configuration the vignette warns against.

## Fixes

* `model_formula()` binds the global environment, so a formula prints without
  the environment tag that varied between builds and leaked into the rendered
  reference pages.

# pilotr 0.2.0

Released 2026-07-15.

## Design specifications

* `pilotr_example()` lists the design specifications shipped with the package, one per
  design family, and returns the path to each for `load_spec()`. They are the same JSON
  files that drive the Python package and the no-code app, and the Python package gains
  the same function.
* The eight per-family specifications now travel inside the package itself, so an
  installed copy can load them without a checkout of the repository.

## Documentation

* A new *Worked examples* article simulates every shipped specification.

# pilotr 0.1.0

Initial release.

## Design specifications and simulation

* `simulate_design()` generates a data set from a portable JSON design specification
  shared with the Python package of the same name, with `build_spec()` composing a
  specification from a flat list of design inputs, `load_spec()` reading one authored
  elsewhere (such as one downloaded from the no-code app) and `spec_json()` serialising
  it back to JSON.
* Response families: Gaussian, lognormal, shifted lognormal, Bernoulli, Poisson, ordinal
  and Beta.
* A shared cross-language random-number generator (`make_rng()`, with the AS 241 inverse
  normal in `as241()`) makes the simulated data bit-identical to the Python package's
  given the same specification and seed.
* A `per_subject` value below 1 or above the number of items is rejected with a clear
  error; the Python package raises on the same inputs.

## Power and precision analysis

* Simulation-based power and design analysis with `power_design()`, `power_mixed()` and
  `power_curve_mixed()`, reporting power alongside the Type S and Type M errors of
  Gelman and Carlin (2014). `power_curve_mixed()` returns an `n_converged` column
  alongside `power` and `type_m`, so the Monte Carlo standard error can be computed over
  the converged replicates.
* Precision and region-of-practical-equivalence (ROPE) design analysis with
  `precision_design()` and `precision_curve()`.
* Every analysis seeds the shared RNG as `seed + (replicate - 1)`, the same indexed-seed
  rule across all power and precision functions, the HPC array script and both languages.
* `power_mixed()` and `power_curve_mixed()` reject a specification without an item unit
  with a clear error instead of returning `NaN` power.
* All simulation-based power and precision analyses take a `workers` argument that
  spreads the Monte Carlo replicates across local cores with base R's `parallel`
  package. Because every replicate seeds the shared RNG from its own index, any worker
  count returns results identical to a serial run, and the sweep functions create their
  worker pool once and reuse it across all sample sizes.

## Modelling bridges and the no-code app

* `model_formula()` and `model_data()` derive the `lmer` formula and modelling data
  frame implied by a specification, and `brms_bridge()` derives a `brms` formula, family
  and priors, mapping every response family to its `brms` counterpart (including
  `Beta()`).
* `generate_r_script()` writes a self-contained, reproducible R script from a
  specification.
* A no-code Shiny application over the same specification, launched with `run_app()`.

## Metadata

* `citation("pilotr")` builds its version note from the package metadata, so it cannot
  drift from `DESCRIPTION`.
