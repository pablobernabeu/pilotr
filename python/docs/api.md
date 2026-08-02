# API reference

The public API, in the groups the [R package's reference
index](https://pablobernabeu.github.io/pilotr/reference/) uses and in the same
order, so a name is found in the same place on either site. Import every name
here from the top-level package, for example `from pilotr import simulate,
power`.

The R package covers more ground: calibration, sweeps, the generated-analysis
emitters and the no-code app have no Python counterpart yet, so the groups they
occupy on that side do not appear below.

## Design specifications

Load and validate the portable JSON specification that everything else reads.
The specification is the same file in either language, and `SPEC_VERSION` is the
version of that format this release understands.

::: pilotr.simulate.load_spec

::: pilotr.validate.validate_spec

::: pilotr.validate.SPEC_VERSION

::: pilotr.examples.pilotr_example

## Simulation

Draw a data set from a specification, returned together with the design
information an analysis needs.

::: pilotr.simulate.simulate

::: pilotr.simulate.Dataset

## Power and design analysis

Estimate power by simulation, at one sample size or across a range of them.

::: pilotr.power.power

::: pilotr.power.power_curve

::: pilotr.power.power_mixed

## Reproducibility

The shared random-number stream and the numerical primitives that make a Python
run and an R run agree bit for bit.

::: pilotr.core.RNG

::: pilotr.core.replicate_seeds

::: pilotr.core.as241

::: pilotr.core.inv_logit
