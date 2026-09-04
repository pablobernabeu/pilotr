## Resubmission

This is a resubmission. The package has not been on CRAN before and the version
is unchanged at 0.3.0.

In response to the comment on writing to the console, `brms_bridge()` no longer
does. It previously wrote the `brms` model it derives to standard output with
`cat()` on every call, whether or not anything was there to read it. It now
returns that model as a `pilotr_bridge` object, carrying the same `formula`,
`family`, `priors` and `code` elements as before, and a registered `print()`
method writes the code and returns the object invisibly. Assigning the result is
therefore silent, as is the one internal use of it inside
`generate_design_analysis()`, which no longer has to divert the output around
itself through a sink; a bare call at the console, or an explicit `print()`,
still shows the model.

Every remaining `cat()` or `print()` call in the package's R code is either
inside one of the two `print()` methods, where console output is the point, or
part of the text of a standalone analysis script that the package emits for the
user to run rather than runs itself. Diagnostics elsewhere use `message()`.

## R CMD check results

Local `R CMD check --as-cran` on a tarball built from the resubmitted sources
(Windows 11 x64, R 4.6.1, pandoc 3.8.3, 2026-08-31):

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
  built tarball (2026-08-31, on the resubmitted sources).
* On GitHub Actions, covering macOS-latest (release), windows-latest (release and
  devel) and ubuntu-latest (release, devel and oldrel-1), each running
  `R CMD check --no-manual --as-cran` (2026-08-20, on the sources as first
  submitted; the change described above is the only one to the R code since).

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
