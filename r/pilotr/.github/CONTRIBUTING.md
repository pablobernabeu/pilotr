<!--
This file is a copy. The original is .github/CONTRIBUTING.md at the repository
root, and any edit belongs there first. pkgdown scans only the package
directory, so without a copy here the R site would have no Community link to
render, and the two have to be kept in step by hand. tests/testthat/test-site.R
checks that they are.
-->

# Contributing to pilotr

Thank you for considering a contribution. Issues and pull requests are welcome,
whether they concern a bug, the documentation or a new feature.

This repository holds pilotr's twin packages, the R package in `r/pilotr/` and
the Python package in `python/`. The no-code app comes in two builds. The Shiny
app in `r/pilotr/inst/app/` is launched by `pilotr::run_app()`. A lighter build
in `app-lite/` runs in the browser without a server and uses the engine files
that `app-lite/engine-files.txt` names. The packages and the builds read the
same portable JSON design specification, documented in `spec/`.

The two packages share a random-number contract. Given the same specification
and seed, they produce bit-identical data wherever they share a maths library,
as they do on Linux. The exception is a rare Beta case under R 4.6.0 or later.
`spec/SPEC.md` sets out that case and what can differ elsewhere. To keep the
contract, a change to the generative core usually has to be made in both
packages.

## Reporting a problem or suggesting a feature

Please open an issue at <https://github.com/pablobernabeu/pilotr/issues>. A
small reproducible example helps a great deal. For a simulation bug, a minimal
design specification, the seed and the output of `sessionInfo()` (R) or your
Python version are ideal.

## Setting up for development

In a clone of the repository, the commands for the R package are:

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
likewise generated. After running both demos, `python python/examples/run_demo.py`
and `Rscript r/pilotr/examples/run_demo.R`, run
`python python/examples/parity_check.py --csv r/pilotr/vignettes/parity-summary.csv`
from the repository root. That script compares CSV files at the precision they
are written with, to a tolerance of 1e-6, so the table summarises the agreement
without proving it.

The R site's favicons in `r/pilotr/pkgdown/favicon/` are built by
`Rscript tools/favicons/build_favicons.R`. When the icon changes, run this
script instead of `pkgdown::build_favicons()`. pkgdown's function would draw
the favicons from the logo, whose wordmark cannot be read at favicon sizes.

Three chunks of the power and precision vignettes fit hundreds of mixed models.
The vignettes therefore show these chunks with `eval = FALSE` and read their
results from three cache files in `r/pilotr/vignettes/`, which hold two tables
and a `power_mixed()` object. Each records the versions of pilotr and lme4 that
produced it. `r/pilotr/tools/vignette-caches.R` writes the caches together with
`cache-fingerprint.txt`, a fingerprint read only by the cache test in
`r/pilotr/tests/testthat/test-vignette-caches.R`. `devtools::test()` runs that
test, and R CMD check skips it. The test fails after a change to the seeds, the
simulated data or the code of those chunks, until the caches are written afresh.

`.github/workflows/refresh-vignette-caches.yaml` writes the caches on GitHub. It
runs when the last commit of a push to a branch other than main has
`[refresh vignette caches]` in its message. Once the workflow is on main, it can
also be started by hand from the Actions tab. The workflow runs the script on
Linux with the current R release, on as many workers as the runner has cores.
When the caches record the version in `r/pilotr/DESCRIPTION` and the cache test
passes, the workflow pushes the files that changed to the branch as a separate
commit and starts CI on that commit. The fingerprint changes only with the
seeds, the data or the chunks, so after a version bump alone that commit holds
three files. Without the workflow, run `Rscript tools/vignette-caches.R [workers]`
from `r/pilotr` and commit the four files it writes. On two workers, a run takes
about 35 minutes.

brms is not a dependency of the package, so the package's test suite cannot
check whether brms accepts the code that `brms_bridge()` and
`generate_design_analysis()` emit. After changing either, run
`Rscript tools/brms/check_bridge.R` with brms installed. It checks every example
specification and a design in each of the two response families no example
uses, and it compiles nothing.

The two builds of the app share the engine in `r/pilotr/R/`, but each has its
own `app.R`. Several rules are written out in both files, so a change to one
`app.R` usually belongs in the other.

The apps' headless tests are scripts outside testthat. Neither
`devtools::test()` nor `devtools::check()` runs them, and of the workflows only
the scheduled dependency-health one does. After changing either app,
`app-lite/engine-files.txt` or an engine file it names, run them from the
repository root with shiny, ggplot2 and callr installed:

```bash
Rscript r/pilotr/tests/test_app_logic.R
Rscript r/pilotr/tests/test_app_server.R
Rscript r/pilotr/tests/test_app_lite_server.R
```

The commands for the Python package are:

```bash
pip install -e "./python[dev]"
cd python
pytest                        # run the test suite
```

The Python documentation site builds with mkdocs from `python/`, with the
plugins and packages that `.github/workflows/site.yml` installs. The guides'
code runs during the build, and the build rewrites `docs/pilotr.bib` whenever
the entry it holds is out of date. Run these commands from the repository root:

```bash
pip install mkdocs-material "mkdocstrings[python]" markdown-exec matplotlib scipy statsmodels pandas
cd python
mkdocs build --strict
```

The harness in `tools/parity/` checks the cross-language contract. Run the
commands below from the repository root. The first two simulate every example
specification and adversarial case from the local R and Python sources. The
third compares the R and Python dumps value for value, against the allowances
in `tools/parity/tolerance.json` and the recorded hashes in
`tools/parity/golden.json`. The next two check that both validators accept and
reject the same specifications with the same messages, and that both curve
solvers agree. The validator check also tests `spec/design.schema.json` against
the validators, for which it needs the jsonschema package. The last three check
the other outputs the two packages share. `json_roundtrip.py` reads back in
Python the numbers R's `spec_json()` writes, `formula_cross.py` compares the
output of the two packages' `model_formula()` and `model_data()`, and
`report_cross.py` compares the text the two `power_report()` functions write.

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

The parity job in `.github/workflows/python-tests.yml` runs all eight. A
failure in its Linux cell fails the workflow. Its Windows cell is
informational, because R and Python on Windows take `exp()`, `log()` and
`pow()` from different maths libraries, and a value that passes through them can
differ in its last bits.

## Conventions

Use British spelling throughout. In Python, the generative core stays
dependency-free, with `scipy`, `statsmodels` and `pandas` as optional extras for
the power tools. The R package keeps its Imports minimal. Changes must keep the
random-number contract and pass the `tools/parity/` harness above. pilotr
handles no credentials and makes no network requests.

## Submitting a pull request

Base your work on `main`, keep the change focused and add or update tests and
documentation alongside the code. The scheduled dependency-health workflow
runs the Python test suite and examples, and the R package's app tests, example
scripts and brms check, against the latest released dependency versions. The
website build exercises the vignettes and guides on every push to `main` that
changes the packages, the app or the specification.

By contributing, you agree that your contribution is licensed under the same MIT
licence as the package, and that you will follow the
[Code of Conduct](CODE_OF_CONDUCT.md).

## Releasing

The R and Python packages are released separately, so their version numbers can
differ. An R release sets `Version` in `r/pilotr/DESCRIPTION`. A Python release
sets the version in `python/pyproject.toml`, `python/pilotr/__init__.py` and
`extra.version` in `python/mkdocs.yml` together. Commit the
`python/docs/pilotr.bib` that a docs build then rewrites. `CITATION.cff` follows
the newer release of the two, with the date it first became public. The README
citations give no version. `python tools/check_versions.py` checks the version
numbers, though not the date, and CI runs it.

On the day CRAN or PyPI accepts a release, push the exact tree submitted, tag
it, push the tag by name (for example `git push origin v0.4.0`, never `--tags`)
and publish the GitHub release, which Zenodo archives.
