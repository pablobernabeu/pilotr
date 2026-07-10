## Submission

This is the first submission of pilotr (version 0.1.0) to CRAN.

## R CMD check results

Local `R CMD check --as-cran` (Windows 11 x64, R 4.6.1, 2026-07-10):

0 errors | 0 warnings | 1 note

The single note is the standard "New submission" note.

## Test environments

* Local: Windows 11 x64, R 4.6.1 (R CMD check --as-cran, 2026-07-10).
* GitHub Actions: windows-latest, macOS-latest, ubuntu-latest
  (release, devel and oldrel-1).

All checks gave the same result: 0 errors, 0 warnings, and only the standard
"New submission" note.

## Notes

* The package contains no compiled code.
* Examples that fit mixed-effects models are wrapped in `\donttest{}` and additionally
  guarded with `requireNamespace("lme4")` / `requireNamespace("lmerTest")`, so they are
  skipped where those suggested packages are unavailable.
* `run_app()` launches an interactive Shiny application; its example is wrapped in
  `\dontrun{}`.
* This is a new package, so there are no reverse dependencies.
