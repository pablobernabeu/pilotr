## Submission

This is the first submission of pilotr to CRAN. The package has not been on CRAN
before, and the version submitted is 0.3.0.

## R CMD check results

Local `R CMD check --as-cran` on a tarball built from the submitted sources
(Windows 11 x64, R 4.6.1, pandoc 3.8.3, 2026-08-21):

0 errors | 0 warnings | 1 note

The note is the new-submission note raised by the CRAN incoming feasibility
check.

    Maintainer: 'Pablo Bernabeu <pcbernabeu@gmail.com>'
    New submission

Examples, examples under `--run-donttest`, the tests, the re-building of the
vignettes and both the PDF and the HTML manual were all checked in that run and
all passed.

## Test environments

Version 0.3.0 was checked in the two environments below.

* Locally on Windows 11 x64 under R 4.6.1, with `R CMD check --as-cran` on the
  built tarball (2026-08-21).
* On GitHub Actions, covering macOS-latest (release), windows-latest (release and
  devel) and ubuntu-latest (release, devel and oldrel-1), each running
  `R CMD check --no-manual --as-cran` (2026-08-20).

All six GitHub Actions runs finished with status OK, meaning no errors, no
warnings and no notes. Those runs disable the CRAN incoming feasibility check,
which is why the new-submission note appears only in the local run.

## Notes

Some aspects of the package are worth flagging for the review.

* The package contains no compiled code.
* Examples that fit mixed-effects models are wrapped in `\donttest{}` and additionally
  guarded with `requireNamespace("lme4")` / `requireNamespace("lmerTest")`, so they are
  skipped where those suggested packages are unavailable.
* `run_app()` launches an interactive Shiny application, so its example is wrapped in
  `\dontrun{}`.
* This is a new package, so there are no reverse dependencies.
