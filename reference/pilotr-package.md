# pilotr: Simulate Experimental and Behavioural Data from a Portable Design Specification

Generative simulation of experimental and behavioural data sets from a
portable JSON design specification shared with the Python package of the
same name. Fixed effect sizes are user-specified, by-subject and by-item
random intercepts and slopes are crossed, and the response families
cover Gaussian, lognormal, shifted lognormal, Bernoulli, Poisson,
ordinal and Beta outcomes. Power and precision-based design analysis run
from the same specification.

## Details

A pilotr workflow begins with a design specification, a plain list
recording the study you plan to run: its groups and conditions, sample
sizes, fixed effect sizes, random-effect standard deviations, and the
response family. Assemble one from a flat list of design inputs with
[`build_spec()`](https://pablobernabeu.github.io/pilotr/reference/build_spec.md),
or read one back from a JSON file with
[`load_spec()`](https://pablobernabeu.github.io/pilotr/reference/load_spec.md).
The package ships one ready-to-run specification per design family, and
[`pilotr_example()`](https://pablobernabeu.github.io/pilotr/reference/pilotr_example.md)
returns their paths.
[`default_response_name()`](https://pablobernabeu.github.io/pilotr/reference/default_response_name.md)
gives the response column that a family uses by default, and
[`spec_json()`](https://pablobernabeu.github.io/pilotr/reference/spec_json.md)
serialises a specification back to JSON for the Python twin or the
no-code app to read.

[`simulate_design()`](https://pablobernabeu.github.io/pilotr/reference/simulate_design.md)
turns a specification into an analysis-ready data frame with one row per
observation. A specification carries its own seed, and both languages
draw from the combined generator built by
[`make_rng()`](https://pablobernabeu.github.io/pilotr/reference/make_rng.md)
on top of the inverse-normal routine
[`as241()`](https://pablobernabeu.github.io/pilotr/reference/as241.md),
so a given specification and seed produce identical data in either
language.

For analysis,
[`model_data()`](https://pablobernabeu.github.io/pilotr/reference/model_data.md)
adds the response column that the model expects, and
[`model_formula()`](https://pablobernabeu.github.io/pilotr/reference/model_formula.md)
derives the maximal mixed-model formula the design implies.
[`brms_bridge()`](https://pablobernabeu.github.io/pilotr/reference/brms_bridge.md)
returns the formula, family and priors for a Bayesian fit.

Design analysis runs from that same specification.
[`power_design()`](https://pablobernabeu.github.io/pilotr/reference/power_design.md)
estimates power for a two-group Gaussian design, together with the Type
S and Type M errors of Gelman and Carlin (2014).
[`power_mixed()`](https://pablobernabeu.github.io/pilotr/reference/power_mixed.md)
does the same for a crossed mixed-effects design, and
[`power_curve_mixed()`](https://pablobernabeu.github.io/pilotr/reference/power_curve_mixed.md)
sweeps sample size to locate where a design becomes adequately powered.
[`precision_design()`](https://pablobernabeu.github.io/pilotr/reference/precision_design.md)
and its curve counterpart
[`precision_curve()`](https://pablobernabeu.github.io/pilotr/reference/precision_curve.md)
report the width of the interval a design buys and the decision
probabilities against a region of practical equivalence.

Two functions round the package off.
[`generate_r_script()`](https://pablobernabeu.github.io/pilotr/reference/generate_r_script.md)
writes a self-contained script that reproduces a simulation, and
[`run_app()`](https://pablobernabeu.github.io/pilotr/reference/run_app.md)
launches the bundled no-code app.

For a worked introduction, see
[`vignette("getting-started", package = "pilotr")`](https://pablobernabeu.github.io/pilotr/articles/getting-started.md).

## See also

Useful links:

- <https://pablobernabeu.github.io/pilotr/>

- <https://github.com/pablobernabeu/pilotr>

- Report bugs at <https://github.com/pablobernabeu/pilotr/issues>

## Author

Pablo Bernabeu, author and maintainer (<pcbernabeu@gmail.com>,
[ORCID](https://orcid.org/0000-0003-1083-2460)).
