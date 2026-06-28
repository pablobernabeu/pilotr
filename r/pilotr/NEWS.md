# pilotr 0.1.0

* Initial release.
* Generative simulation of experimental and behavioural data from a portable JSON design
  specification shared with the Python package of the same name, using a shared
  cross-language random-number generator for bit-identical output.
* Response families: Gaussian, shifted lognormal, Bernoulli, Poisson, ordinal and Beta.
* Simulation-based power and precision/ROPE design analysis, reporting power alongside the
  Type S and Type M errors of Gelman and Carlin (2014).
* A no-code Shiny application over the same specification, launched with `run_app()`.
