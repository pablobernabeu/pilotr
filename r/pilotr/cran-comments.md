## Submission

This release fixes the test failure reported for 0.3.0 on the no-long-double
check.

    -- Failure ('test-core.R:296:3'): spec_json round-trips a coefficient exactly --
    Expected `as.numeric(back$fixed$coefficients$cond)` to be identical to `1/3`.
    Differences:
      `actual`: 0.333333333333333259
    `expected`: 0.333333333333333315

The fault was in the package rather than in the test, and it was a real one:
the saved design specification no longer held the coefficient the user had set.
`spec_json()` wrote the document by asking `jsonlite` for 17 significant digits
and then shortening each number in the resulting text, which meant reading the
number back with `as.numeric()`. Without long doubles `as.numeric()` accumulates
the mantissa in a double before applying the decimal exponent, and 17
significant digits overflow the 53-bit mantissa on the way, so
`"0.33333333333333331"` reads back one unit in the last place below `1/3`. The
shortened form then recorded that wrong value.

Each number is now formatted from the double itself, before `jsonlite` sees the
document, so nothing is read back from text the package has just written. A
coefficient of `1/3` is written at 16 significant digits rather than 17, which
is both the shorter form and the one such a build reads back exactly. The same
helper feeds `generate_r_script()`, where it could return nothing at all for a
finite value and so embed `Inf` in place of a coefficient. It now falls back to
17 significant digits instead.

I do not have a no-long-double build to hand, so the fix was checked by
reproducing in R the accumulation `R_strtod()` performs when `LDOUBLE` is
`double`. Under that arithmetic the 17-digit form gives 0.333333333333333259,
the value reported above, and the 16-digit form the package now writes gives
`1/3` exactly. Two regression tests were added: one reads every number
`spec_json()` writes under that same arithmetic, the other covers a label
written to imitate the marker the new code uses internally.

## R CMD check results

Local `R CMD check --as-cran` on a tarball built from the submitted sources
(Windows 11 x64, R 4.6.1, pandoc 3.10, 2026-09-17):

0 errors | 0 warnings | 1 note

The note comes from the CRAN incoming feasibility check and is about the
interval since the last release.

    Maintainer: 'Pablo Bernabeu <pcbernabeu@gmail.com>'

    Days since last update: 5

0.3.0 was published five days ago and this release exists only to correct the
fault that the no-long-double check found in it, so the short interval is the
reason for the submission rather than an oversight.

Examples, examples under `--run-donttest`, the tests, the re-building of the
vignettes and both the PDF and the HTML manual were all checked in that run and
all passed.

## Test environments

Version 0.3.1 was checked locally on Windows 11 x64 under R 4.6.1, with
`R CMD check --as-cran` on the built tarball (2026-09-17).

The repository's GitHub Actions workflow covers macOS-latest (release),
windows-latest (release and devel) and ubuntu-latest (release, devel and
oldrel-1), each running `R CMD check --no-manual --as-cran`, and runs on the
commit carrying this change.

## Notes

* The package contains no compiled code.
* Examples that fit mixed-effects models are wrapped in `\donttest{}` and additionally
  guarded with `requireNamespace("lme4")` / `requireNamespace("lmerTest")`, so they are
  skipped where those suggested packages are unavailable.
* `run_app()` launches an interactive Shiny application, so its example is wrapped in
  `\dontrun{}`.
* There are no reverse dependencies.
