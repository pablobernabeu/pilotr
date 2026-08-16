# Changelog

All notable changes to this project are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed
- **A Poisson mean beyond the sampler's reach was returned as the iteration cap.** The
  inverse-CDF walk starts from `exp(-mean)`, which underflows to exactly zero once the
  mean passes about 746 (a poisson intercept of 7 already implies a mean of exp(7),
  about 1097), after which every simulated count came back as 1000000 while reporting
  success. Both engines now refuse such a mean, naming `exp(eta)` and the offending
  value, with message text character-for-character identical across the twins; a linear
  predictor large enough to overflow `math.exp` itself folds into the same refusal
  instead of surfacing as `OverflowError`. Feasible means are untouched and the parity
  dumps are unchanged.
- **A true effect of exactly zero broke both engines, and the package recommends that
  input.** `design_conditions()` deliberately produces a null condition so a run can show
  how often it declares something when there is nothing to find; feeding it through
  `power()` raised `ZeroDivisionError` in Python and returned `Inf` in R, with Type S
  quietly degenerating to "the estimate is positive". Type S and Type M are now
  undefined when the true effect is zero or unknown, in both engines — the guard the
  mixed-effects path already carried, applied to the other three sites.
- **Four small twin divergences closed**: a non-whole seed now truncates identically in
  both engines, a whole `spec_version` is read the same however it was written, a
  non-object unit is reported rather than crashing R with a base error, and
  `replicate_seeds()` is annotated as returning integers. Message text is now
  character-for-character identical across the two engines in each case.

### Changed
- **The Bayesian design-analysis record carries the rule that produced it** — the decision
  thresholds, the convergence gate, a fingerprint of the specification and the pilotr
  version — and the aggregator refuses to pool replicates that disagree on any of them.
  Two array runs with different ROPEs previously aggregated into one table with nothing
  to tell them apart.
- **The parity harness anchors only the IEEE-754-exact cases.** The golden anchor's
  criterion was "zero ulp allowance", which admitted cases whose bytes pass through
  `exp()`, `log()` or `pow()` whenever the platforms measured so far happened to agree
  on them. `tolerance.json` now classifies those cases as transcendental explicitly:
  they stay gated by the cross-language comparison, at zero ulp unless an allowance is
  recorded, and `golden.json` pins only the Gaussian cases. The validator cross-check
  (`tools/parity/validate_cross.py`) also runs in CI now, and its default `Rscript` is
  the one on `PATH` rather than a path on the author's machine.

### Added
- **A test that the packaged example specifications match the repository's own.**
  `spec/examples/*.json` is canonical and both packages carry a mirror so an
  installed copy can reach it, but nothing enforced the mirror. The existing test
  could not: a stale packaged copy still loads and simulates perfectly well, it
  simply describes a different design from the one the repository documents. The
  new test compares the bytes, and is twinned with `test-examples.R`.
- **CI tests what users install, not only the checkout.** The matrix gains Python
  3.14. A new job installs the built wheel into a bare environment outside the
  repository and calls `pilotr_example()` there, then resolves and parses one
  specification: that is the only check that can catch the example mirror going
  missing from the distribution, since an editable install finds the repository
  copy regardless. A second new job installs the declared minimum dependency
  versions of the extras. A weekly schedule runs the suite when nobody has
  pushed; the workflow's paths filters mean a quiet month otherwise sees no run
  at all.

<!-- The date is when the 0.3.0 work landed on main, not a release date: 0.3.0 has
     not been tagged yet. Set it to the tag date when the release is cut, along with
     date-released in CITATION.cff. -->
## [0.3.0] - 2026-08-01

Two of the changes below alter numbers that earlier versions produced. Both are
deliberate corrections, and neither can be made without moving the output, so
install 0.2.1 to reproduce output from 0.2.1.

### Added

- `validate_spec()` checks a design specification against the portable schema and
  against the cross-field rules the schema cannot express, and `load_spec()` and
  `simulate()` now call it by default. A strict draft-07 schema had shipped since
  0.1 with no code path consulting it, and several ways of getting a
  specification wrong produced plausible data rather than an error. A mistyped
  coefficient key resolved to no column and so silently set that effect to zero,
  which generates exactly the data of a null design and reports success. A
  response parameter left over from another family was ignored. Both are now
  refused. Both functions take a `validate` argument, so a replicate loop can
  validate once and then skip the check.
- `SPEC_VERSION`, and the `spec_version` field it negotiates. A specification
  with no such field is a 0.2 specification. One that uses a 0.3 feature has to
  declare 0.3, because a 0.2 implementation reads it differently and generates
  different data without complaint, and one declaring a version newer than this
  implementation understands is refused rather than partly read.
- `replicate_seeds()` is exported, so a hand-written loop or a cluster array task
  can seed its replicates by the same rule the power functions use.
- The `exgaussian` response family, the registered model family for
  reaction-time work, taking `sigma` and `beta`. It is written in brms's
  parameterisation, in which `mu` is the mean, so a specification and the model
  fitted to it agree on what the intercept means.
- A predictor can now vary by `"observation"`, drawing a fresh value for every
  row, for a quantity that varies trial by trial.
- A predictor can be drawn from a uniform distribution through `dist` with `min`
  and `max`. A uniform costs the same single draw as a normal, so switching
  between them does not move the random stream.
- A predictor can carry a `reliability`, which simulates imperfect measurement.
  The latent value drives the linear predictor and any random slope keyed on the
  predictor, while the contaminated observed value is what appears in the
  returned data, which is what an analyst would actually have measured.

### Changed

- The power functions no longer seed replicate `i` with `seed + (i - 1)`.
  Consecutive seeds are not independent streams in this generator, and the
  measured consequence was severe: the first draw of replicate `i` correlated
  0.95 with the first draw of replicate `i + 1`. Seeds now come from the shared
  generator through `replicate_seeds()`, which brings that correlation to -0.02.
  Every number `power()`, `power_curve()` and `power_mixed()` produce is
  therefore different from 0.2.1, and the first replicate no longer uses the
  specification's own seed.
- A random slope keyed on an interaction, such as `"z_cosine:z_ISI"`, now reaches
  the data. It was accepted, sized into the covariance and drawn from the random
  stream, and then looked up as though it were a plain column, which returned
  nothing, so it was silently discarded. Any design that used one was
  anti-conservative and its data have changed.
- The linear predictor is accumulated as an explicit left fold in the order set
  out in `spec/SPEC.md`. Floating-point addition is not associative, and R folded
  the terms one at a time while Python summed them and added the total, so the
  two disagreed on 63.8% of rows in a matched reproduction over 200,000 rows.
- No built-in summation is used anywhere in the generative core. CPython's `sum()`
  has applied Neumaier compensation to floats since 3.12 and base R's `sum()`
  accumulates in 80-bit long double. Both are more accurate than a plain double
  fold, in different ways, so an inner product of three terms or more could land
  on different doubles in the two languages. Every inner product is now an
  explicit fold.
- The cube in the Gamma sampler is written as two multiplications rather than
  `** 3`. Python calls the library `pow()` while R special-cases small integer
  exponents, and over 200,000 draws in the sampler's range the two disagreed on a
  third of inputs by up to 6 ulp. That value decides a rejection step, so it also
  changed how many draws were consumed, which accounted for every difference in
  the `beta_proportion` example.
- The random-effect covariance is bracketed as `(sd_i · sd_j) · r_ij`, matching
  R. Multiplication is commutative but not associative in floating point, and the
  difference propagated through the Cholesky factor into every random effect
  drawn.

### Fixed

- A random-effect covariance that is not positive definite is now an error naming
  the grouping factor and the column at which the factorisation failed. The
  failing pivot was previously clamped at zero and the result returned in
  silence, which produced random effects whose standard deviations were several
  times those requested.
- A correlation naming a random-effect term that does not exist is reported
  against the entry the user wrote, rather than failing as a lookup error deep
  inside the covariance code.
- `varies_by` is validated against the three unit names that exist. Anything
  other than `"subject"` was previously read as item-level, so a predictor
  declared to vary by `"trial"` was silently given one value per item.

## [0.2.1] - 2026-07-23

### Changed

- The guides load a bundled design specification through `pilotr_example()`
  rather than by a path relative to the repository. The path resolved only
  inside a checkout, so a reader who had installed pilotr and copied the line
  met a missing-file error.

## [0.2.0] - 2026-07-15

### Added

- `pilotr_example()` lists the design specifications shipped with the package,
  one per design family, and returns the path to each for `load_spec()`. The
  eight JSON files now travel inside the wheel, so the *Worked examples* page and
  any installed copy can reach them. This matches the R twin's `pilotr_example()`.

## [0.1.0] - 2026-07-10

The first public version of the Python package. It simulates experimental and
behavioural data from a portable JSON design specification, producing output that
is bit-identical to the R package of the same name given the same specification
and seed.

### Added

- `simulate` generates a dataset from a design specification supplied as a plain
  dictionary or a JSON file, with `load_spec` reading a specification authored
  elsewhere, such as one downloaded from the no-code app. A `per_subject` value
  below 1 or above the number of items is rejected with a clear `ValueError`; the
  R package raises on the same inputs.
- Response families for Gaussian outcomes, reaction times, accuracy, counts,
  ordinal responses and Beta-distributed proportions.
- Power and design analysis through `power`, `power_curve` and `power_mixed`,
  with `scipy` and `statsmodels` supplied as optional extras so the generative
  core stays dependency-free. A specification without an item unit is rejected
  by `power_mixed` with a clear `ValueError` rather than yielding `nan` power.
- A cross-language random-number contract that reproduces the R package's data to
  full floating-point precision. Every Monte Carlo replicate seeds the shared RNG
  as `seed + (replicate - 1)`, the same indexed-seed rule as the R package across
  all analyses.
- `power`, `power_mixed` and `power_curve` take a `workers` argument that spreads
  the Monte Carlo replicates across local processes with
  `concurrent.futures.ProcessPoolExecutor`. Because every replicate seeds the
  shared RNG from its own index, any worker count returns results identical to a
  serial run, and `power_curve` starts one process pool and reuses it across the
  whole sweep.
- A documentation site whose guides execute their examples at build time, so the
  tables and figures shown are real `pilotr` output.

[Unreleased]: https://github.com/pablobernabeu/pilotr/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/pablobernabeu/pilotr/releases/tag/v0.3.0
[0.2.1]: https://github.com/pablobernabeu/pilotr/releases/tag/v0.2.1
[0.2.0]: https://github.com/pablobernabeu/pilotr/releases/tag/v0.2.0
[0.1.0]: https://github.com/pablobernabeu/pilotr/releases/tag/v0.1.0
