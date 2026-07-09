# pilotr 0.2.0

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
