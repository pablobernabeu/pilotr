## Submission

This is the first submission of pilotr (version 0.2.1) to CRAN.

## R CMD check results

Local `R CMD check --as-cran` (Windows 11 x64, R 4.6.1, 2026-07-23):

0 errors | 0 warnings | 2 notes

The single note is the standard "New submission" note.

The local machine has no pandoc on the check subprocess PATH, so the local run
also reports that 'README.md' and 'NEWS.md' cannot be checked. That note is an
artefact of this machine and does not arise where pandoc is present, including
the GitHub Actions runs above.


## Test environments

* Local: Windows 11 x64, R 4.6.1 (R CMD check --as-cran, 2026-07-23).
* GitHub Actions: windows-latest (release and devel), macOS-latest (release),
  ubuntu-latest (release, devel and oldrel-1).

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
