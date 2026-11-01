# Changelog

All notable changes to this project are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.4.0] - YYYY-MM-DD

### Upgrading from 0.3.0

#### Changed results

- `power_mixed` fits every coefficient and only the declared random intercepts
  and slopes, each as an independent variance component, where it fitted one
  contrast with by-subject and by-item intercepts and slopes whatever was
  declared. Declared correlations draw a warning.
- Each point of `power_curve` draws its own replicates, so the results at every
  point change. Shared seeds made `solve_curve` intervals too narrow.
  `seeds="shared"` restores the 0.3.0 seeding.
- `power_mixed` fits by Powell's method, with L-BFGS as a fallback, where
  statsmodels' default optimisers often stopped short of the REML optimum and
  inflated standard errors. Power usually rises and, in small designs, can
  exceed the power the R package gets from `lmer()` and Satterthwaite tests.
- `power_mixed` analyses a `lognormal` response on the log scale, as it did a
  `shifted_lognormal` one, where it fitted the raw values and so reported
  estimates in the response's units.
- For a `bernoulli`, `poisson`, `ordinal` or `beta` response, `power_mixed`
  warns and returns the mean estimate and Type M as `nan`, where it compared a
  response-scale estimate with a link-scale true value.
- `n_converged` from `power_mixed` counts only the fits that are neither
  singular nor flagged with a warning, where it counted every fit that returned,
  as `n_returned` now does.
- Power and the mean estimate from `power_mixed` leave out fits without a finite
  estimate and p-value, which 0.3.0 counted as not significant and averaged.
- With `validate=False`, `simulate` allocates two or more factors between one
  unit jointly, one block per combination of levels, where each factor took its
  own blocks, which coincided or overlapped and left cells empty.
- Rounded responses and Poisson counts are computed as in R, where CPython's
  `round()` and the bracketing of each term of the Poisson recurrence could give
  another value at a decimal tie or, rarely, another count. The shipped examples
  simulate the same data.
- A Beta response that came out as exactly 0 or 1, which brms's `Beta()`
  refuses, becomes the nearest double inside (0, 1) or, for a rounded response,
  the nearest value inside it at its number of decimals. No other value changes.
- A Beta response whose Gamma shape is exactly 1/2 is drawn as the R twin draws
  it, where `pow()` gave another double for some uniforms on Linux. The last
  bits of a few rows of such a design change.

#### Changed interface

- `power` and `power_curve` refuse a design whose rows are correlated by an item
  unit, a within factor or a grouping factor besides `subject`, with
  `NotImplementedError`, where `power` t-tested every row as independent and
  overstated power.
- `power_mixed` returns each focal effect's results under `effects`. Reading one
  of the top-level keys from `power` to `type_m_mcse`, which hold the results of
  the effect 0.3.0 tested, issues a `FutureWarning` (see Deprecated). `backend`
  names the test and the direction of a one-sided test.
- `validate_spec` refuses two or more factors between the same unit unless the
  specification declares `"spec_version": "0.4"`, where such factors fell into
  overlapping blocks that left cells empty. One between factor whose levels are
  the cells works in every version.
- `validate_spec` refuses colliding names, such as a factor called `subject` or
  a level listed twice, which moved an effect or overwrote a column, and blank
  names. It also refuses a `random.item` entry in a design without items, which
  was dropped.
- `validate_spec` refuses a factor that sets both `vary_within` and `between`,
  which gave each unit one row per level, all at one level. It also refuses a
  correlation that pairs a term with itself, which rescaled its variance, and
  two correlations for one pair.
- `power_mixed` refuses a design that declares no random effects with
  `NotImplementedError`, where, with an item unit, it fitted by-subject and
  by-item components whatever the specification declared and, without one,
  raised `ValueError`.
- `power` refuses a design with fewer than two subjects at a level of the
  between factor, and `power_curve` a grid size below 4, where three subjects in
  all left the pooled variance one degree of freedom and two gave a power of 0.
- Results gain keys. `power` returns `alternative`, and `power_mixed` returns
  the fit counts listed under Added. The power functions also return the Monte
  Carlo errors listed under Added, and each `power_curve` record also gains
  `type_s`.
- In a data frame built from `power_curve` records, `type_m`, `n_sims` and
  `n_significant` sit further right, so code that selects them by position gets
  other columns.
- `simulate` refuses a Poisson mean above about 708, where counts were biased or
  intermittently refused from a mean of about 738 upwards.
- `simulate` refuses extreme linear predictors with `ValueError` and the same
  message as the R twin, where it raised a bare `OverflowError` for an
  overflowing lognormal or shifted lognormal response and `ZeroDivisionError`
  for a Beta response.
- `validate_spec` raises `ValueError` for every malformed shape, as documented,
  where some raised `TypeError` or `AttributeError`.
- `validate_spec` refuses a seed beyond ±(2^53 − 1), which R and Python can read
  as different numbers, and a `round` of 0 for a Beta response, which rounded
  every response to 0 or 1.
- `load_spec` refuses a key repeated within one object, which R and Python read
  differently, even with `validate=False`, and so does every function given a
  path.
- `load_spec` reads a specification as UTF-8 whatever the locale, where, on
  Windows, it used the locale encoding, which misread or refused UTF-8 labels
  outside ASCII. A file in another encoding with such labels usually raises
  `UnicodeDecodeError`.
- `validate_spec` refuses a `random` of `0` or `false`, which was read as no
  random effects, and `"item": null` among the units, which passed and then
  raised `TypeError` in `simulate`.
- `SPEC_VERSION` is `"0.4"`, where it was `"0.3"`.
- `pilotr.core.beta_draw(rng, a, b)`, which the documented API never included,
  is replaced by `beta_response(rng, eta, phi)`.

#### Changed output

- `Dataset.to_csv` writes UTF-8 whatever the locale, where, on Windows, it
  wrote the locale encoding, which raised `UnicodeEncodeError` for some labels
  outside ASCII and wrote others in bytes that UTF-8 readers refuse.
- `Dataset.to_csv` quotes a field holding a comma, a double quote or a newline,
  where a comma shifted the later fields of its row and a double quote or a
  newline could break the row.
- A NumPy number in a specification reaches the simulated data as a Python
  number, where, under NumPy 2, `to_csv` wrote it as `np.float64(...)`.
- A within factor whose `vary_within` omits a unit of the design warns (see
  Deprecated).
- `validate_spec` reads a unit given as `[]` as an empty object, as R does, and
  so reports its missing `n`, where it said that the unit must be an object.

### Added

- Specification version 0.4 allocates several between factors on one unit
  jointly, one block per combination of levels with the first-listed factor
  varying slowest. pilotr 0.3.0 refuses a specification that declares it, and
  earlier releases do not check the version.
- `counterbalance`, a specification-0.4 field on a within factor, gives each
  subject one condition per item, rotating the lists as in a Latin square, and
  adds a `list` column to the simulated data. `power_mixed` analyses such a
  design.
- Every rate carries its Monte Carlo error: the standard error and 95% Wilson
  interval of power (`power_mcse`, `power_lo`, `power_hi`), the standard error
  of Type M (`type_m_mcse`) and a Wilson interval for Type S (`type_s_lo`,
  `type_s_hi`). `power`, `power_mixed` and `power_curve` report them all.
- `power`, `power_mixed` and `power_curve` gain `alternative`, which plans a
  one-sided test with `"greater"` or `"less"`. The default, `"two.sided"`, is
  the test of earlier releases, and scipy's spelling `"two-sided"` is refused.
- `power_mixed` analyses any design that declares a random effect, with or
  without items, with any number of within factors and with an extra grouping
  factor such as a site.
- `model_formula`, a port of the R function, writes the model a specification
  implies in lme4's notation, maximal, zero-correlation or
  random-intercepts-only as `structure` says, and `model_data`, also a port,
  builds the columns that model reads.
- `power_mixed` gains `focal`, which chooses the effects to test and their true
  values, by default every coefficient with the value the specification gives
  it. A contrast without a coefficient is left out of the model, so a null
  effect is planned with a coefficient of 0. For a `bernoulli`, `poisson`,
  `ordinal` or `beta` response, the Type S of an interaction is `nan`.
- `power_mixed` reports `n_attempted`, `n_returned`, `n_singular` and
  `n_warning`, as the R twin does, judging a fit singular by lme4's test.
- `power_report` writes a sample-size justification in Markdown from a power
  result and, optionally, a solved curve, with the items of Lakens (2022,
  <https://doi.org/10.1525/collabra.33267>) and the simulation accounting of
  Morris, White and Crowther (2019, <https://doi.org/10.1002/sim.8086>). The R
  twin writes the same text from the same counts.

### Changed

- `power`, `power_curve` and `power_mixed` validate the specification once
  before their replicates, where each replicate validated it again.

### Deprecated

- A within factor whose `vary_within` omits a unit of the design warns, since
  pilotr crosses a within factor with every unit whatever the list says. List
  every unit instead, or make the factor between `item` if items carry it. From
  specification version 0.4, such a list is refused.
- The top-level keys of a `power_mixed` result from `power` to `type_m_mcse`
  will be removed in the next release. Read them from `effects` instead.

### Fixed

- `load_spec` reads a specification as R's jsonlite does, so a one-element array
  where a single value belongs, `[]` for an empty `coefficients`, `slopes` or
  `correlations` and a byte-order mark are accepted, where each raised an error
  or was refused.
- A `null` in an optional field other than a unit or an entry of `random`, such
  as `predictors`, `correlations` or `spec_version`, reads as an absent field,
  in a file or a dict, where it raised an error or was refused.
- `simulate` accepts, as R does, a whole number written as a float, such as
  `"n": 64.0`, a single ordinal threshold written as a bare number and, with
  `validate=False`, a `round` that is not a whole number, where each raised
  `TypeError`.
- `power` reads a missing coefficient for the tested contrast as a true effect
  of 0, with Type S and Type M undefined, where it raised `KeyError`.
- `power_mixed` returns no estimate for a replicate whose response has no
  logarithm, such as a shifted lognormal value rounded down to the shift, where
  the whole call stopped.
- NumPy numbers and `pathlib` paths are accepted wherever a number or a path is,
  and pandas data frames wherever a curve is, where a NumPy integer was refused
  as not a whole number and a path as not a JSON object. `validate_spec` accepts
  a path, `workers` any whole number and `solve_curve` a NumPy `n` as a
  replicate count.

### Documentation

- The cross-language guarantee is stated as it holds: R and Python simulate
  bit-identical data where they share a maths library, as on Linux, except a
  Beta response whose Gamma shape is exactly 1/3 or 1/4 under R 4.6.0 or later.
  It was claimed for the Gaussian family on every platform and for any rounded
  response. The documentation also says that one language's data can differ
  across operating systems in the same way, where it said that a specification
  and seed always give the same data.
- The specification sets out how units are allocated to levels and clusters,
  with worked encodings for a 2 × 2 between design, a two-list counterbalanced
  design and randomisation within clusters, and gives every sampling rule in
  pseudocode, correcting its account of the draws the Gamma sampler consumes.
  Its JSON schema requires each family's parameters and refuses unknown
  top-level fields and unused response parameters. It accepts the single-value
  forms of `vary_within` and `thresholds`, as strict validation does.
- The `solve_curve` docstring says that when an effect varies from item to item,
  power at a fixed number of items levels off below 1, however many subjects are
  added, so a curve over subjects can be solved far from the true sample size.
- The `solve_curve` docstring reports the coverage of its interval over 1,000
  curves a setting, where it reported coverage in 35 of 36 checks on curves
  whose points shared their seeds. It says that such a curve, including one
  built from runs that each start from the same seed, gets too narrow an
  interval, and that a finer grid helps only while its points are independent.
- The `inv_logit` docstring gives its range as [0, 1], since it returns exactly
  1 from an argument of about 36.7 and exactly 0 below about −745.1, where it
  gave (0, 1).
- The power guide gains sections on one-sided tests, the model and test each
  family gets and a sample-size justification with `power_report`.
- The power guide says that its crossed mixed-effects example is fully crossed
  and that a counterbalanced study has half the observations and is planned with
  `counterbalance` or the two-list encoding. The specification says that
  simulated data are complete, so a study that expects to lose a proportion p
  of its participants recruits N / (1 − p) to analyse N.
- The documentation states which random effects `power_mixed` fits, each as an
  independent variance component, and which features only the R package
  offers, where this package's README and home page called it the R package's
  feature-parity twin and the project README put the two languages at
  capability parity.
- The install commands in the README and Get started carry no trailing comments
  and so paste into zsh, where pip read the comments of the 0.3.0 commands as
  requirements and stopped with 'Invalid requirement'.

## [0.3.0] - 2026-08-27

Two of the changes below alter numbers that earlier versions produced. Both are
deliberate corrections, and neither can be made without moving the output, so
install 0.2.1 to reproduce output from 0.2.1.

### Added

- `validate_spec()` checks a design specification against the portable schema and
  against the cross-field rules the schema cannot express, and `load_spec()` and
  `simulate()` now call it by default. A strict draft-07 schema had shipped since
  0.1 with no code path consulting it, and several ways of getting a
  specification wrong produced plausible data and no error at all. A mistyped
  coefficient key resolved to no column and so silently set that effect to zero,
  which generates exactly the data of a null design and reports success. A
  response parameter left over from another family was ignored. Both are now
  refused. Both functions take a `validate` argument, so a replicate loop can
  validate once and then skip the check.
- `SPEC_VERSION`, and the `spec_version` field it negotiates. A specification
  with no such field is a 0.2 specification. One that uses a 0.3 feature has to
  declare 0.3, because a 0.2 implementation reads it differently and generates
  different data without complaint, and one declaring a version newer than this
  implementation understands is refused outright.
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

- `solve_curve` and `target_n` solve a simulated design curve for the value that meets a
  target. The package computed the whole curve and then handed back the last and most
  consequential step, the number that goes into a preregistration, at replicate counts
  where neighbouring points are not significantly different. The documentation drew a dashed
  line at 0.80 and left the reader to judge the crossing by eye. `target_n` now takes the
  records `power_curve` already returns and reports the sample size at which power reaches
  the target, with a confidence interval on it, rounded up to whole subjects. `solve_curve`
  is the general form and reads a curve on any axis through the same column names.

  The fit is a binomial probit regression of the decision rate against the swept value,
  weighted by the replicates behind each point and inverted by the delta method that R's
  `MASS::dose.p` applies to a fitted `glm`. The probit follows from the design itself:
  under the normal approximation to a two-group comparison the probit of power is linear
  in the square root of the sample size. Checked against R's `stats::power.t.test` over 36
  combinations of effect size, target power and grid shape, at 400 replicates a point, the
  solved sample size sat within 2.9% of the analytic answer on average and 6.9% at worst,
  against 3.4% and 8.9% for a logit fitted the same way.
  `tools/calibration/solve_curve_calibration.R` runs that comparison and writes every
  figure quoted for it to a file beside itself. Where the two-parameter model does not
  describe the curve, the interval is widened by the heterogeneity factor of probit
  analysis, reported as `dispersion`.

  Nothing extrapolates. A curve that does not reach the target within the sizes it swept is
  refused, with the range it did cover, and so is a fit that solves past the end of the
  sweep, a curve with no trend to invert, and a slope that cannot be told from zero. Every
  refusal message is character-for-character identical with the R twin's, checked by
  `tools/parity/solve_cross.py`.

- `power_curve` reports the replicate counts behind each point. Its records gain
  `n_sims` and `n_significant` alongside `n_subject`, `power` and `type_m`. A power estimate
  over 1000 replicates should weigh more in a solve than one over 50, and the count is what
  says which it is. Without it a curve could not be handed to `solve_curve` at all.

- The packaged example specifications are now tested against the repository's own copies.
  `spec/examples/*.json` is canonical and both packages carry a mirror so an
  installed copy can reach it, but nothing enforced the mirror. The existing test
  could not: a stale packaged copy still loads and simulates perfectly well, it
  simply describes a different design from the one the repository documents. The
  new test compares the bytes, and is twinned with `test-examples.R`.
- CI now tests the built distribution as well as the checkout. The matrix gains
  Python 3.14. A new job installs the built wheel into a bare environment outside the
  repository and calls `pilotr_example()` there, then resolves and parses one
  specification: that is the only check that can catch the example mirror going
  missing from the distribution, since an editable install finds the repository
  copy regardless. A second new job installs the declared minimum dependency
  versions of the extras. A weekly schedule runs the suite when nobody has
  pushed, since the workflow's paths filters mean a quiet month otherwise sees no
  run at all.

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

- The Bayesian design-analysis record carries the rule that produced it, namely the
  decision thresholds, the convergence gate, a fingerprint of the specification and the
  pilotr version, and the aggregator refuses to pool replicates that disagree on any of
  them. Two array runs with different ROPEs previously aggregated into one table with
  nothing to tell them apart.
- The parity harness anchors only the IEEE-754-exact cases. The golden anchor's
  criterion was "zero ulp allowance", which admitted cases whose bytes pass through
  `exp()`, `log()` or `pow()` whenever the platforms measured so far happened to agree
  on them. `tolerance.json` now classifies those cases as transcendental explicitly:
  they stay gated by the cross-language comparison, at zero ulp unless an allowance is
  recorded, and `golden.json` pins only the Gaussian cases. The validator cross-check
  (`tools/parity/validate_cross.py`) also runs in CI now, and its default `Rscript` is
  whichever one is on `PATH`, where it used to be a path on the author's machine.
- The parity contract now has a third gate, at a relative tolerance where the other two
  work in ulps. `solve_curve` is an iterative fit, and it calls `exp` and the normal
  distribution function at every step, so a golden hash would pin one platform's maths
  library exactly as `tolerance.json` says such cases must not be. It is kept out of
  `golden.json` and out of the dump harness, and gated by `tools/parity/solve_cross.py`,
  which puts fixed curves through both engines and compares the solved value, its
  interval and the text of every refusal. The allowance is 1e-9 relative, against a largest
  measured disagreement of 5.0e-15.

### Fixed

- A random-effect covariance that is not positive definite is now an error naming
  the grouping factor and the column at which the factorisation failed. The
  failing pivot was previously clamped at zero and the result returned in
  silence, which produced random effects whose standard deviations were several
  times those requested.
- A correlation naming a random-effect term that does not exist is reported
  against the entry the user wrote, where it used to fail as a lookup error deep
  inside the covariance code.
- `varies_by` is validated against the three unit names that exist. Anything
  other than `"subject"` was previously read as item-level, so a predictor
  declared to vary by `"trial"` was silently given one value per item.

- A Poisson mean beyond the sampler's reach was returned as the iteration cap. The
  inverse-CDF walk starts from `exp(-mean)`, which underflows to exactly zero once the
  mean passes about 746 (a poisson intercept of 7 already implies a mean of exp(7),
  about 1097), after which every simulated count came back as 1000000 while reporting
  success. Both engines now refuse such a mean, naming `exp(eta)` and the offending
  value, with message text character-for-character identical across the twins. A linear
  predictor large enough to overflow `math.exp` itself folds into the same refusal, so it
  never reaches the caller as an `OverflowError`. Feasible means are untouched and the
  parity dumps are unchanged.
- A true effect of exactly zero broke both engines, and the package recommends that
  input. `design_conditions()` deliberately produces a null condition so a run can show
  how often it declares something when there is nothing to find. Feeding it through
  `power()` raised `ZeroDivisionError` in Python and returned `Inf` in R, with Type S
  quietly degenerating to "the estimate is positive". Type S and Type M are now
  undefined when the true effect is zero or unknown, in both engines, which is the guard
  the mixed-effects path already carried, applied to the other three sites.
- Four small twin divergences are closed. A non-whole seed now truncates identically in
  both engines, a whole `spec_version` is read the same however it was written, a
  non-object unit is now reported as such, where R used to crash with a base error, and
  `replicate_seeds()` is annotated as returning integers. Message text is
  character-for-character identical across the two engines in each case.

## [0.2.1] - 2026-07-23

### Changed

- The guides load a bundled design specification through `pilotr_example()`. The
  path relative to the repository that they used before resolved only inside a
  checkout, so a reader who had installed pilotr and copied the line met a
  missing-file error.

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
  below 1 or above the number of items is rejected with a clear `ValueError`, and
  the R package raises on the same inputs.
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

[Unreleased]: https://github.com/pablobernabeu/pilotr/compare/v0.4.0...HEAD
[0.4.0]: https://github.com/pablobernabeu/pilotr/releases/tag/v0.4.0
[0.3.0]: https://github.com/pablobernabeu/pilotr/releases/tag/v0.3.0
[0.2.1]: https://github.com/pablobernabeu/pilotr/releases/tag/v0.2.1
[0.2.0]: https://github.com/pablobernabeu/pilotr/releases/tag/v0.2.0
[0.1.0]: https://github.com/pablobernabeu/pilotr/releases/tag/v0.1.0
