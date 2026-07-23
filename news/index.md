# Changelog

## pilotr 0.2.1

Released 2026-07-23.

### Documentation

- A package-level help topic,
  [`?pilotr`](https://pablobernabeu.github.io/pilotr/reference/pilotr-package.md),
  introduces the workflow and groups every export by the stage it
  belongs to.
- The reference examples that cut the replicate count so they run
  quickly now say so, and each `n_sims` argument documents the Monte
  Carlo error and a usable recommendation for real planning.
- The
  [`precision_curve()`](https://pablobernabeu.github.io/pilotr/reference/precision_curve.md)
  example passes a region of practical equivalence clearly inside the
  effect it declares. It previously fell back to a default equal to that
  effect, the configuration the vignette warns against.

### Fixes

- [`model_formula()`](https://pablobernabeu.github.io/pilotr/reference/model_formula.md)
  binds the global environment, so a formula prints without the
  environment tag that varied between builds and leaked into the
  rendered reference pages.

## pilotr 0.2.0

Released 2026-07-15.

### Design specifications

- [`pilotr_example()`](https://pablobernabeu.github.io/pilotr/reference/pilotr_example.md)
  lists the design specifications shipped with the package, one per
  design family, and returns the path to each for
  [`load_spec()`](https://pablobernabeu.github.io/pilotr/reference/load_spec.md).
  They are the same JSON files that drive the Python package and the
  no-code app, and the Python package gains the same function.
- The eight per-family specifications now travel inside the package
  itself, so an installed copy can load them without a checkout of the
  repository.

### Documentation

- A new *Worked examples* article simulates every shipped specification.

## pilotr 0.1.0

Initial release.

### Design specifications and simulation

- [`simulate_design()`](https://pablobernabeu.github.io/pilotr/reference/simulate_design.md)
  generates a data set from a portable JSON design specification shared
  with the Python package of the same name, with
  [`build_spec()`](https://pablobernabeu.github.io/pilotr/reference/build_spec.md)
  composing a specification from a flat list of design inputs,
  [`load_spec()`](https://pablobernabeu.github.io/pilotr/reference/load_spec.md)
  reading one authored elsewhere (such as one downloaded from the
  no-code app) and
  [`spec_json()`](https://pablobernabeu.github.io/pilotr/reference/spec_json.md)
  serialising it back to JSON.
- Response families: Gaussian, lognormal, shifted lognormal, Bernoulli,
  Poisson, ordinal and Beta.
- A shared cross-language random-number generator
  ([`make_rng()`](https://pablobernabeu.github.io/pilotr/reference/make_rng.md),
  with the AS 241 inverse normal in
  [`as241()`](https://pablobernabeu.github.io/pilotr/reference/as241.md))
  makes the simulated data bit-identical to the Python package’s given
  the same specification and seed.
- A `per_subject` value below 1 or above the number of items is rejected
  with a clear error; the Python package raises on the same inputs.

### Power and precision analysis

- Simulation-based power and design analysis with
  [`power_design()`](https://pablobernabeu.github.io/pilotr/reference/power_design.md),
  [`power_mixed()`](https://pablobernabeu.github.io/pilotr/reference/power_mixed.md)
  and
  [`power_curve_mixed()`](https://pablobernabeu.github.io/pilotr/reference/power_curve_mixed.md),
  reporting power alongside the Type S and Type M errors of Gelman and
  Carlin (2014).
  [`power_curve_mixed()`](https://pablobernabeu.github.io/pilotr/reference/power_curve_mixed.md)
  returns an `n_converged` column alongside `power` and `type_m`, so the
  Monte Carlo standard error can be computed over the converged
  replicates.
- Precision and region-of-practical-equivalence (ROPE) design analysis
  with
  [`precision_design()`](https://pablobernabeu.github.io/pilotr/reference/precision_design.md)
  and
  [`precision_curve()`](https://pablobernabeu.github.io/pilotr/reference/precision_curve.md).
- Every analysis seeds the shared RNG as `seed + (replicate - 1)`, the
  same indexed-seed rule across all power and precision functions, the
  HPC array script and both languages.
- [`power_mixed()`](https://pablobernabeu.github.io/pilotr/reference/power_mixed.md)
  and
  [`power_curve_mixed()`](https://pablobernabeu.github.io/pilotr/reference/power_curve_mixed.md)
  reject a specification without an item unit with a clear error instead
  of returning `NaN` power.
- All simulation-based power and precision analyses take a `workers`
  argument that spreads the Monte Carlo replicates across local cores
  with base R’s `parallel` package. Because every replicate seeds the
  shared RNG from its own index, any worker count returns results
  identical to a serial run, and the sweep functions create their worker
  pool once and reuse it across all sample sizes.

### Modelling bridges and the no-code app

- [`model_formula()`](https://pablobernabeu.github.io/pilotr/reference/model_formula.md)
  and
  [`model_data()`](https://pablobernabeu.github.io/pilotr/reference/model_data.md)
  derive the `lmer` formula and modelling data frame implied by a
  specification, and
  [`brms_bridge()`](https://pablobernabeu.github.io/pilotr/reference/brms_bridge.md)
  derives a `brms` formula, family and priors, mapping every response
  family to its `brms` counterpart (including `Beta()`).
- [`generate_r_script()`](https://pablobernabeu.github.io/pilotr/reference/generate_r_script.md)
  writes a self-contained, reproducible R script from a specification.
- A no-code Shiny application over the same specification, launched with
  [`run_app()`](https://pablobernabeu.github.io/pilotr/reference/run_app.md).

### Metadata

- `citation("pilotr")` builds its version note from the package
  metadata, so it cannot drift from `DESCRIPTION`.
