# About

## Citing pilotr

If you use pilotr in published work, please cite it:

> Bernabeu, P. (2026). *pilotr: Simulate experimental and behavioural data from a portable
> design specification* (Version 0.1.0) [Computer software].
> https://doi.org/10.5281/zenodo.21266313

```bibtex
@software{bernabeu2026pilotr,
  author  = {Pablo Bernabeu},
  title   = {{pilotr}: Simulate experimental and behavioural data from a portable design specification},
  year    = {2026},
  version = {0.1.0},
  doi     = {10.5281/zenodo.21266313},
  url     = {https://doi.org/10.5281/zenodo.21266313}
}
```

[Download the BibTeX entry](pilotr.bib){ download="pilotr.bib" }

In R, `citation("pilotr")` returns the same reference.

## The developer

pilotr is developed by [Pablo Bernabeu](https://pablobernabeu.github.io), a researcher in
the Department of Education at the University of Oxford. His work spans cognitive psychology,
neuroscience, linguistics, education and research methods, with hands-on experience of
behavioural experiments, EEG, corpus analysis, computational modelling and statistics. He
develops open, reproducible research software in R and Python, and is a Fellow of the Software
Sustainability Institute. pilotr and its
[R twin](https://pablobernabeu.github.io/pilotr/) share one design specification, keeping a
simulation reproducible across both languages. His
[ORCID record](https://orcid.org/0000-0003-1083-2460) lists his other work.

## Licence

pilotr is released under the MIT licence, reproduced in full on the
[licence page](licence.md). The licence covers the Python and R packages, the no-code app
and the design specification alike.

## Versioning and archival

Each release is tagged on GitHub and archived on Zenodo. The concept DOI,
[10.5281/zenodo.21266313](https://doi.org/10.5281/zenodo.21266313), always resolves to the
latest archived version, so a citation stays current without naming a version. The
[changelog](changelog.md) records what changed in each release.

## Contributing and support

Bugs and feature requests are best raised on the
[GitHub issues page](https://github.com/pablobernabeu/pilotr/issues). The
[contributing guide](https://github.com/pablobernabeu/pilotr/blob/main/.github/CONTRIBUTING.md)
describes the development setup and the conventions the repository follows, including the
cross-language random-number contract that any change to the generative core must keep
intact.
