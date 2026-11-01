<!--
This file is a copy. The original is .github/CONTRIBUTING.md at the repository
root, and any edit belongs there first. pkgdown scans only the package
directory, so without a copy here the R site would have no Community link to
render, and the two have to be kept in step by hand. tests/testthat/test-site.R
checks that they are.
-->

# Contributing to pilotr

Thank you for considering a contribution. Issues and pull requests are both
welcome, whether they fix a bug, improve the documentation or add a feature.

pilotr ships as twin packages, R and Python, from this one repository. The R
package lives in `r/pilotr/` and the Python package in `python/`. The no-code
app comes in two builds: the Shiny app in `r/pilotr/inst/app/`, which
`pilotr::run_app()` launches, and a lighter build in `app-lite/` that runs in
the browser without a server, on the engine files `app-lite/engine-files.txt`
names. All of them consume the same portable JSON design
specification (documented in `spec/`). The two packages share a random-number
contract so that, given the same specification and seed, they produce
bit-identical data wherever they share a maths library, as on Linux, outside a
rare Beta case under R 4.6.0 or later. `spec/SPEC.md` sets out that case and
what can differ elsewhere. A change to the generative core therefore usually
needs to land in both packages, keeping that contract intact.

## Reporting a problem or suggesting a feature

Please open an issue at <https://github.com/pablobernabeu/pilotr/issues>. A
small reproducible example helps a great deal; for a simulation bug, a minimal
design specification plus the seed is ideal, along with the output of
`sessionInfo()` (R) or your Python version.

## Setting up for development

Working from a clone of the repository, the R package:

```r
pak::pak(c("devtools", "cran::roxygen2@8.0.0", "testthat", "spelling"))
devtools::document("r/pilotr")   # regenerate man/ and NAMESPACE after editing roxygen
devtools::test("r/pilotr")       # run the test suite
devtools::check("r/pilotr")      # a full R CMD check
devtools::build_readme("r/pilotr")  # regenerate README.md after editing README.Rmd
```

Install the roxygen2 release that `Config/roxygen2/version` in
`r/pilotr/DESCRIPTION` records, which is the one CI pins. A newer release
rewrites that field when `devtools::document()` runs, and the CI job that checks
`NAMESPACE` then fails.

`devtools::check()` also runs a spelling test, as CI does, and stops at any word
found neither in the British English dictionary nor in `r/pilotr/inst/WORDLIST`.
`spelling::spell_check_package("r/pilotr")` lists those words. Correct a
misspelling where it occurs, and add a name or a term of art to the word list.
CRAN and a plain `R CMD check` skip the test.

`r/pilotr/README.md` is generated from `README.Rmd` and is also the website's home
page, so edit the `.Rmd` and knit it rather than editing the Markdown directly.
The cross-language parity table in `r/pilotr/vignettes/parity-summary.csv` is
likewise generated, from the repository root, by
`python python/examples/parity_check.py --csv r/pilotr/vignettes/parity-summary.csv`,
which needs both demos, `python python/examples/run_demo.py` and
`Rscript r/pilotr/examples/run_demo.R`, to have been run first. That table
compares CSV files at the precision
they are written with, to a tolerance of 1e-6, so it summarises the agreement
without proving it. The R site's favicons in `r/pilotr/pkgdown/favicon/`
are built by `Rscript tools/favicons/build_favicons.R`. Run it, not
`pkgdown::build_favicons()`, when the icon changes, since pkgdown's function
would draw them from the logo, whose wordmark cannot be read at favicon sizes.

brms is not a dependency of the package, so its test suite cannot ask brms
whether it accepts the code that `brms_bridge()` and `generate_design_analysis()`
emit. After changing either, run `Rscript tools/brms/check_bridge.R` with brms
installed. It checks every shipped example and a design in each of the two
response families no example uses, and it compiles nothing.

The two builds of the app share the engine in `r/pilotr/R/`, but each has its
own `app.R`, and several rules are written out in both, so a change to one
`app.R` usually belongs in the other. Their headless tests are scripts outside
testthat. Neither `devtools::test()` nor `devtools::check()` runs them, and of
the workflows only the scheduled dependency-health one does. After changing
either app, `app-lite/engine-files.txt` or an engine file it names, run them from
the repository root with shiny, ggplot2 and callr installed:

```bash
Rscript r/pilotr/tests/test_app_logic.R
Rscript r/pilotr/tests/test_app_server.R
Rscript r/pilotr/tests/test_app_lite_server.R
```

The Python package:

```bash
pip install -e "./python[dev]"
cd python
pytest                        # run the test suite
```

The Python documentation site builds with mkdocs from `python/`, with the
plugins and packages that `.github/workflows/site.yml` installs. The guides'
code runs during the build, and the build rewrites `docs/pilotr.bib` whenever
the entry it holds is out of date. From the repository root:

```bash
pip install mkdocs-material "mkdocstrings[python]" markdown-exec matplotlib scipy statsmodels pandas
cd python
mkdocs build --strict
```

The cross-language contract is checked by the harness in `tools/parity/`. From
the repository root, the first two commands below simulate every shipped and
adversarial case from the local R and Python sources. The third compares the
two dumps value for value, against the allowances in
`tools/parity/tolerance.json` and the recorded hashes in
`tools/parity/golden.json`. The next two check that both validators accept and
refuse the same specifications with the same messages, and that both curve
solvers agree. The validator check also holds `spec/design.schema.json` to the
validators, for which it needs the jsonschema package. The last three hold the
other outputs the twins share: `json_roundtrip.py` reads back in Python the
numbers R's `spec_json()` writes, `formula_cross.py` compares the two
`model_formula()` and `model_data()`, and `report_cross.py` compares the text
the two `power_report()` functions write.

```bash
Rscript tools/parity/run_r.R
python tools/parity/run_py.py
python tools/parity/compare.py
python tools/parity/validate_cross.py
python tools/parity/solve_cross.py
python tools/parity/json_roundtrip.py
python tools/parity/formula_cross.py
python tools/parity/report_cross.py
```

The parity job in `.github/workflows/python-tests.yml` runs all eight. Its
Linux cell gates the workflow. Its Windows
cell is informational, because there R and Python take `exp()`, `log()` and
`pow()` from different maths libraries, and a value that passes through them can
differ in its last bits.

## Conventions

British spelling throughout. The generative core stays dependency-free in
Python (`scipy`, `statsmodels` and `pandas` are optional extras for the power
tools), and the R package keeps its Imports minimal. Given the same
specification and seed, the two languages produce bit-identical data wherever
they share a maths library, and the `tools/parity/` harness above must stay
green. pilotr handles no credentials and makes no network requests.

## Submitting a pull request

Base your work on `main`, keep the change focused, and add or update tests and
documentation alongside the code. A scheduled continuous-integration workflow
runs the Python test suite and examples, and the R package's app tests, example
scripts and brms check, against the latest released dependency versions, and the
website build exercises the vignettes and guides on every push to `main` that
changes the packages, the app or the specification.

By contributing you agree that your contribution is licensed under the same MIT
licence as the package, and that you will follow the
[Code of Conduct](CODE_OF_CONDUCT.md).

## Releasing

The R and Python packages are released separately, so their version numbers can
differ. An R release sets `Version` in `r/pilotr/DESCRIPTION`. A Python release
sets the version in `python/pyproject.toml`, `python/pilotr/__init__.py` and
`extra.version` in `python/mkdocs.yml` together, after which a docs build
rewrites `python/docs/pilotr.bib` for committing. `CITATION.cff` follows the
newer release of the two, with the date it first became public. The README
citations give no version. `python tools/check_versions.py` checks the version
numbers, though not the date, and CI runs it.

On the day CRAN or PyPI accepts a release, push the exact tree submitted, tag it
by name (for example `git push origin v0.4.0`, never `--tags`) and publish the
GitHub release, which Zenodo archives.
