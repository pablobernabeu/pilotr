# Contributing to pilotr

Thank you for considering a contribution. Issues and pull requests are both
welcome, whether they fix a bug, improve the documentation or add a feature.

pilotr ships as twin packages, R and Python, from this one repository. The R
package lives in `r/pilotr/` and the Python package in `python/`, with the
no-code app in `app-lite/`. All three consume the same portable JSON design
specification (documented in `spec/`), and the two packages share a
random-number contract so that, given the same specification and seed, they
produce bit-identical data. A change to the generative core therefore usually
needs to land in both packages, keeping that contract intact.

## Reporting a problem or suggesting a feature

Please open an issue at <https://github.com/pablobernabeu/pilotr/issues>. A
small reproducible example helps a great deal; for a simulation bug, a minimal
design specification plus the seed is ideal, along with the output of
`sessionInfo()` (R) or your Python version.

## Setting up for development

Working from a clone of the repository, the R package:

```r
pak::pak(c("devtools", "roxygen2", "testthat", "spelling"))
devtools::document("r/pilotr")   # regenerate man/ and NAMESPACE after editing roxygen
devtools::test("r/pilotr")       # run the test suite
devtools::check("r/pilotr")      # a full R CMD check
```

The Python package:

```bash
pip install -e "./python[dev]"
cd python
pytest                        # run the test suite
python examples/parity_check.py   # confirm cross-language agreement
```

## Conventions

British spelling throughout. The generative core stays dependency-free in
Python (`scipy`, `statsmodels` and `pandas` are optional extras for the power
tools), and the R package keeps its Imports minimal. Outputs are bit-identical
across the two languages given the same specification and seed;
`python/examples/parity_check.py` exercises that contract and must stay green.
pilotr handles no credentials and makes no network requests.

## Submitting a pull request

Base your work on `main`, keep the change focused, and add or update tests and
documentation alongside the code. Continuous integration runs both test suites
on every push.

By contributing you agree that your contribution is licensed under the same MIT
licence as the package, and that you will follow the
[Code of Conduct](CODE_OF_CONDUCT.md).
