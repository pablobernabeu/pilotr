# pilotr 0.2.0

* Precision replicates (`precision_design()`, `precision_curve()` and the HPC array script)
  seed the shared RNG as `seed + (replicate - 1)`, the same indexed-seed rule as the power
  functions; they previously used `seed + replicate`. This pre-release alignment keeps the
  seed rule identical across all analyses and languages, and it changes the results of
  precision runs performed before the alignment.
* `power_curve_mixed()` returns an `n_converged` column alongside `power` and `type_m`, so
  the Monte Carlo standard error can be computed over the converged replicates.
* `power_mixed()` and `power_curve_mixed()` reject a specification without an item unit with
  a clear error instead of returning `NaN` power.
* `simulate_design()` rejects a `per_subject` value below 1 or above the number of items with
  a clear error; the Python package raises on the same inputs.
* `build_spec()` supports the lognormal family: `sigma` is carried through and the default
  response name matches the shifted lognormal's.
* `brms_bridge()` maps the beta family to `brms`'s `Beta()`.
* `citation("pilotr")` builds its version note from the package metadata, so it can no longer
  drift from `DESCRIPTION`.
* All simulation-based power and precision analyses (`power_design()`, `power_mixed()`,
  `power_curve_mixed()`, `precision_design()` and `precision_curve()`) gain a `workers`
  argument that spreads the Monte Carlo replicates across local cores with base R's
  `parallel` package. Because every replicate seeds the shared cross-language RNG from its
  own index, any worker count returns results identical to a serial run.
* The sweep functions create their worker pool once and reuse it across all sample sizes,
  so a high-resolution curve pays the cluster start-up cost only once.

# pilotr 0.1.0

* Initial release.
* Generative simulation of experimental and behavioural data from a portable JSON design
  specification shared with the Python package of the same name, using a shared
  cross-language random-number generator for bit-identical output.
* Response families: Gaussian, shifted lognormal, Bernoulli, Poisson, ordinal and Beta.
* Simulation-based power and precision/ROPE design analysis, reporting power alongside the
  Type S and Type M errors of Gelman and Carlin (2014).
* A no-code Shiny application over the same specification, launched with `run_app()`.
