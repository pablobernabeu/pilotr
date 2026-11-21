# pilotr (Python)

<!-- badges: start -->
[![PyPI](https://img.shields.io/pypi/v/pilotr.svg)](https://pypi.org/project/pilotr/)
[![python-tests](https://github.com/pablobernabeu/pilotr/actions/workflows/python-tests.yml/badge.svg)](https://github.com/pablobernabeu/pilotr/actions/workflows/python-tests.yml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/license/MIT)
[![DOI 10.5281/zenodo.21266313](https://zenodo.org/badge/DOI/10.5281/zenodo.21266313.svg)](https://doi.org/10.5281/zenodo.21266313)
<!-- badges: end -->

The Python implementation of pilotr simulates experimental and behavioural data from a
portable JSON design specification.

This is the twin of [the R package](https://pablobernabeu.github.io/pilotr/r/) of the same name
and contains its generative core. Both packages use the same design specification and
random-number generator, so the same specification and seed produce identical data in either
language wherever the two share a maths library, as on Linux, except in a rare Beta case under
R 4.6.0 or later that the
[cross-language guide](https://pablobernabeu.github.io/pilotr/python/cross-language/) describes.
Elsewhere, a value that passes through `exp()`, `log()` or `pow()` can differ in its last bits,
as the same guide sets out.

This package provides simulation, validation, the analysis model a specification implies,
simulation-based power, the sample-size report and the curve solvers. Precision/ROPE design
analysis, power curves for mixed-effects designs, calibration and sweeps, the `lme4` reference
backend for mixed-effects power, the `brms` bridge, the analysis-script emitters, the helpers
that build specifications, and the app launcher are in the R package alone.

## Install

pilotr requires Python 3.9 or later. The first command below installs the generative core, which
is pure Python and has no dependencies. The second adds scipy, for `power` and `power_curve`, and
the third adds scipy, statsmodels and pandas, for `power_mixed` and `model_data`. The double
quotes keep zsh, the default shell on macOS, from reading the brackets as a filename pattern.

```bash
pip install pilotr
pip install "pilotr[power]"
pip install "pilotr[mixed]"
```

For development, install from a checkout of the repository instead. Clone
<https://github.com/pablobernabeu/pilotr> and run `pip install .` (or `pip install ".[dev]"`)
from the `python` directory.

## Quick start

A design specification is an ordinary dictionary, so the shortest example needs nothing beyond
the package itself.

```python
from pilotr import simulate

spec = {
    "name": "two_group", "seed": 2024,
    "units": {"subject": {"n": 64}},
    "factors": [{"name": "group", "levels": ["control", "treatment"],
                 "contrasts": {"effect": [-0.5, 0.5]}, "between": "subject"}],
    "fixed": {"intercept": 100, "coefficients": {"effect": 5}},
    "response": {"family": "gaussian", "name": "score", "sigma": 10},
}

data = simulate(spec)               # 64 rows
```

Given the same specification and seed, this reproduces the data the R package generates, to the
last bit wherever the two share a maths library.

## Try it without installing

A serverless build of the [no-code app](https://pablobernabeu.github.io/pilotr/app/) runs in
your browser, with no installation required and no data uploaded.

## Learn more

The [Python documentation site](https://pablobernabeu.github.io/pilotr/python/) walks through the
quick start, the response families, power and design analysis, the specification format and the
API reference. The [project README](https://github.com/pablobernabeu/pilotr) describes the
toolkit as a whole, including the R package, the no-code app and the specification format.

## Citation

If pilotr contributes to published work, please cite it, with the version of the package you
used (`pilotr.__version__`). The R and Python packages are released separately, so their version
numbers can differ.

> Bernabeu, P. (2026). *pilotr: Simulate experimental and behavioural data from a portable
> design specification*. <https://doi.org/10.5281/zenodo.21266313>

The [About page](https://pablobernabeu.github.io/pilotr/python/about/) gives the same citation,
with the Python package's version and a BibTeX entry, and a short note on the developer. The
repository also includes
[`CITATION.cff`](https://github.com/pablobernabeu/pilotr/blob/main/CITATION.cff), which GitHub's
'Cite this repository' button reads.

## Licence

MIT. See [LICENSE](https://github.com/pablobernabeu/pilotr/blob/main/python/LICENSE).

## Contributing

Issues and pull requests are welcome. The
[contributing guide](https://github.com/pablobernabeu/pilotr/blob/main/.github/CONTRIBUTING.md)
describes the development setup and the conventions the package follows, and everyone taking
part is asked to honour the
[Code of Conduct](https://github.com/pablobernabeu/pilotr/blob/main/.github/CODE_OF_CONDUCT.md).
