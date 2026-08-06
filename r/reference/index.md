# Package index

## Design specifications

Write, validate and load the portable JSON specification that everything
else reads.

- [`build_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/build_spec.md)
  : Build a design specification from a flat list of design inputs
- [`validate_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/validate_spec.md)
  : Validate a design specification
- [`load_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/load_spec.md)
  : Load a design specification from a JSON file
- [`spec_json()`](https://pablobernabeu.github.io/pilotr/r/reference/spec_json.md)
  : Serialise a design specification to pretty-printed JSON
- [`spec_from_model()`](https://pablobernabeu.github.io/pilotr/r/reference/spec_from_model.md)
  : Build a design specification from a fitted mixed model
- [`pilotr_example()`](https://pablobernabeu.github.io/pilotr/r/reference/pilotr_example.md)
  : Example design specifications shipped with pilotr

## Simulation

Draw a data set from a specification, together with the model that
specification implies.

- [`simulate_design()`](https://pablobernabeu.github.io/pilotr/r/reference/simulate_design.md)
  : Simulate a data set from a design specification
- [`model_data()`](https://pablobernabeu.github.io/pilotr/r/reference/model_data.md)
  : Build the modelling data frame from a simulated data set and its
  specification
- [`model_formula()`](https://pablobernabeu.github.io/pilotr/r/reference/model_formula.md)
  : Derive the lmer formula implied by a specification
- [`default_response_name()`](https://pablobernabeu.github.io/pilotr/r/reference/default_response_name.md)
  : Default response-column name for a family

## Power and design analysis

Estimate power by simulation, at one sample size or across a range of
them.

- [`power_design()`](https://pablobernabeu.github.io/pilotr/r/reference/power_design.md)
  : Simulation-based power and design analysis for a two-group Gaussian
  design
- [`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md)
  : Simulation-based power and design analysis for a mixed-effects
  design
- [`power_curve_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_curve_mixed.md)
  : Power curve over sample size for a mixed-effects design
- [`print(`*`<pilotr_power>`*`)`](https://pablobernabeu.github.io/pilotr/r/reference/print.pilotr_power.md)
  : Print a simulation-based power result

## Precision and equivalence

Size a study for the width of an interval rather than for a significance
test.

- [`precision_design()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md)
  : Precision and ROPE design analysis at a fixed sample size
- [`precision_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_curve.md)
  : Precision and ROPE curve over sample size

## Calibration and sweeps

Tune a specification towards a target, and vary one of its fields across
a range.

- [`calibrate_response()`](https://pablobernabeu.github.io/pilotr/r/reference/calibrate_response.md)
  : Rescale a design to a target total variance
- [`response_variance()`](https://pablobernabeu.github.io/pilotr/r/reference/response_variance.md)
  : Variance components of the linear predictor
- [`design_conditions()`](https://pablobernabeu.github.io/pilotr/r/reference/design_conditions.md)
  : Build a grid of fixed-effect coefficient sets
- [`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md)
  : Sweep an analysis over one field of a design specification

## Generated analysis code

Emit runnable, self-contained code for the analysis a specification
implies.

- [`generate_r_script()`](https://pablobernabeu.github.io/pilotr/r/reference/generate_r_script.md)
  : Generate a self-contained, reproducible R script from a
  specification
- [`generate_design_analysis()`](https://pablobernabeu.github.io/pilotr/r/reference/generate_design_analysis.md)
  : Generate a Bayesian design-analysis script from a specification
- [`brms_bridge()`](https://pablobernabeu.github.io/pilotr/r/reference/brms_bridge.md)
  : Derive a brms formula, family, and priors from a design spec

## Reproducibility

The shared random-number stream that makes an R run and a Python run
agree bit for bit.

- [`make_rng()`](https://pablobernabeu.github.io/pilotr/r/reference/make_rng.md)
  : Create a shared cross-language random-number generator
- [`replicate_seeds()`](https://pablobernabeu.github.io/pilotr/r/reference/replicate_seeds.md)
  : Seeds for the replicates of a Monte Carlo run
- [`as241()`](https://pablobernabeu.github.io/pilotr/r/reference/as241.md)
  : Inverse normal cumulative distribution function (Wichura's AS 241)

## App

The code-free interface to the same workflow.

- [`run_app()`](https://pablobernabeu.github.io/pilotr/r/reference/run_app.md)
  : Launch the pilotr no-code app

## Package overview

What the package does, and where to start.

- [`pilotr`](https://pablobernabeu.github.io/pilotr/r/reference/pilotr-package.md)
  [`pilotr-package`](https://pablobernabeu.github.io/pilotr/r/reference/pilotr-package.md)
  : pilotr: Simulate Experimental and Behavioural Data from a Portable
  Design Specification
