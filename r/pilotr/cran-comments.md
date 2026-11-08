## Submission

This is an update from pilotr 0.3.1, which CRAN published on 2026-09-20, to
0.4.0.

0.4.0 changes results and interfaces. Some functions return different values
for the same input: `power_mixed()`, `precision_design()` and their curves fit
'bernoulli' and 'poisson' responses with `lme4::glmer()`, and each point of a
sweep draws its own replicates, so power, interval widths and Type M change.
Some specifications that 0.3.1 accepted are refused, results gain elements and
columns, and printed output changes. NEWS.md lists each of these changes under
'Upgrading from the previous release (0.3.1)', in three parts: 'Changed
results', 'Changed interface' and 'Changed output'.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

GitHub Actions, run 42427181969 on 2026-11-07, for commit d6541f0 on the
release branch. Pull request #24 started the run, so it checked out
0fe482d, the merge of d6541f0 into main at 0221556. That merge differs from
d6541f0 only in .github/workflows/publish.yml, outside the package. The three
files of cached mixed-model results in vignettes/ are the ones that the
refresh-vignette-caches workflow regenerated and committed as cde86bd, so that
they record 0.4.0. None of the commits since cde86bd changed them, the chunks
that produced them or the functions those chunks call. The package sources of
this submission are those of d6541f0. Only this file, which `.Rbuildignore`
keeps out of the tarball, has changed since. Each leg built the tarball with
`--no-manual` and `--compact-vignettes=gs+qpdf` and ran
`R CMD check --as-cran --no-manual` on it, without the CRAN incoming checks,
which r-lib/actions turns off. Each leg ended with 'Status: OK':

* ubuntu-latest (Ubuntu 24.04.5 LTS), R 4.6.2
* ubuntu-latest (Ubuntu 24.04.5 LTS), R-devel (2026-11-06 r90845)
* ubuntu-latest (Ubuntu 24.04.5 LTS), R 4.5.3
* windows-latest (Windows Server 2022 x64, build 26100), R 4.6.2
* windows-latest (Windows Server 2022 x64, build 26100), R-devel (2026-11-06
  r90845)
* macos-latest (macOS Tahoe 26.7.2, aarch64-apple-darwin23), R 4.6.2

The test suite passed in every leg, with no failures or warnings. 34 tests skip
everywhere because they read repository files that R CMD check cannot reach,
such as the shared specification schema, the browser app and the vignette
sources, and one because pkgdown is not installed. Under R 4.5.3, one more
skips because it checks how R 4.6.0 and later compute integer powers.

The same run installed the Suggests floors that DESCRIPTION declares (ggplot2
3.4.0 and testthat 3.1.5) on R oldrel-3 under ubuntu-latest. There the check,
run without the vignettes, ended with no errors or warnings, and its two notes
concern packages that the floors job does not install.

## Notes

* The package contains no compiled code.
* lme4 and lmerTest are in Suggests. The mixed-model functions stop with a
  message naming them when they are missing, and the examples that fit mixed
  models are wrapped in `\donttest{}` and guarded with `requireNamespace()`, so
  they are skipped where those packages are unavailable.
* shiny, ggplot2, future and promises are used by the no-code app that
  `run_app()` launches, and callr by the app's Verify button, which asks for
  callr to be installed when it is absent. The app is interactive, so the
  example of `run_app()` is wrapped in `\dontrun{}`.
* MASS, codetools and spelling are used by the tests only. Each test that needs
  one of them, or callr, skips when it is absent. testthat runs the tests, and
  knitr and rmarkdown build the vignettes.
* Three vignette chunks that fit hundreds of mixed models are shown with
  `eval = FALSE`, and the vignettes read their results from files in
  vignettes/, so that the vignettes build in seconds.

## Downstream dependencies

There are no reverse dependencies on CRAN.
