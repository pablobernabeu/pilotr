# pilotr <span class="mrd-lang">(Python)</span>

<p class="mrd-tagline">Simulate experimental and behavioural data from a portable JSON design specification.</p>

This is the feature-parity twin of [the R package](https://pablobernabeu.github.io/pilotr/r/) of
the same name. The two share the design specification and the random-number generator, so the
same specification and seed produce identical data in either language, bit for bit apart from a
documented tolerance of a few units in the last place where an unrounded response family applies
`exp()` or `log()` to the linear predictor.

[Get started](getting-started.md){ .md-button .md-button--primary }
[Try the no-code app](https://pablobernabeu.github.io/pilotr/app/){ .md-button }

A design specification names the units, the factors and their contrasts, the fixed effect
sizes, any random effects and a response family. From that one file pilotr draws a data set,
estimates power and precision by simulation, then emits the analysis code the design implies.
The same file drives the R package, the browser-based no-code app and this Python package, so a
design built by pointing and clicking can be run unchanged in a script, and a design written in
a script can be handed to a collaborator who has installed nothing.

The generative core is pure Python and dependency-free. `scipy`, `statsmodels` and `pandas` are
optional extras, needed only by the power tools. The [Get started](getting-started.md) guide
covers installing pilotr and simulating a first data set.

## Where to go next

The guides follow the workflow, from choosing a response family to sizing a study and checking
that the two languages agree.

- [Response families](families.md): Gaussian, lognormal, reaction times, accuracy, counts,
  ordinal, Beta and the ex-Gaussian.
- [Worked examples](examples.md): one ready-to-run design per family.
- [Power and design analysis](power.md): `power`, `power_curve` and `power_mixed`.
- [Cross-language reproducibility](cross-language.md): how the same spec gives identical data
  in R and Python.
- [Specification format](specification.md): the JSON spec and the cross-language RNG contract.
- [API reference](api.md): every public function and class.

Archived on Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21266313.svg)](https://doi.org/10.5281/zenodo.21266313)
