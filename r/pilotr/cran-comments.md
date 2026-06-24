## Submission

This is the first submission of pilotr to CRAN.

## R CMD check results

Local `R CMD check --as-cran` (Windows 11 x64, R 4.6.0):

0 errors | 0 warnings | 1 note

The single note is the standard "New submission" note.

## Test environments

* Local: Windows 11 x64, R 4.6.0, with pandoc 3.8.3.

Before acceptance I will also run win-builder (R-release and R-devel) and a macOS
environment (mac-builder / R-hub).

## Notes

* The package contains no compiled code.
* Examples that fit mixed-effects models are wrapped in `\donttest{}` and additionally
  guarded with `requireNamespace("lme4")` / `requireNamespace("lmerTest")`, so they are
  skipped where those suggested packages are unavailable.
* `run_app()` launches an interactive Shiny application; its example is wrapped in
  `\dontrun{}`.
* This is a new package, so there are no reverse dependencies.
