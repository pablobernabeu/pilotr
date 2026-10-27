# Changelog

All notable changes to this project are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Specification version 0.4 allocates several between factors on one unit jointly, so a 2 × 2
  between-subjects design over 40 subjects has four cells of 10. The units fall into one block
  per combination of levels, with the first-listed factor varying slowest, and nothing is drawn.
  A specification using this must declare `"spec_version": "0.4"`, a version pilotr 0.3.0
  refuses to read. One that does not still gets the refusal described under Changed, which names
  both remedies. A 0.4 specification whose `vary_within` omits a unit of the design is refused.
  `SPEC_VERSION` is now "0.4", the newest version this release reads. The R package's
  `build_spec()` and `spec_from_model()` declare the lowest version a specification needs, so
  the specifications the no-code app exports stay readable by pilotr 0.3.0 unless they use a 0.4
  feature.
- `counterbalance` on a within factor gives each subject one condition per item, rotating the
  lists across subjects and items, as in a Latin-square design. Planning such a study with a
  fully crossed specification doubled the observations per participant and overstated power.
  Subject s meets item t in condition ((s − 1) + (t − 1)) mod C, where C counts the conditions,
  or the combinations of levels when several factors are counterbalanced, and nothing is drawn.
  The field belongs to specification version 0.4, needs an item unit and is refused with
  `per_subject`. A warning, in the R twin's words, says when the lists are unequal, because the
  subjects or the items are not a multiple of C. The condition stays a within factor, so
  `power_mixed` analyses a counterbalanced design, as it does the two-list encoding, which stays
  documented for releases that read no version above 0.3.
- The simulated data of a counterbalanced specification name each subject's list in a column
  `list`, placed after `subject` and `item`. Subject s follows list ((s − 1) mod C) + 1, and
  nothing is drawn for it, so a specification without a counterbalanced factor simulates as
  before. Such a design keeps the name `list` for that column, and `validate_spec` refuses a
  factor, predictor, grouping factor, contrast column or response of that name in it, in the R
  twin's words.
- Every rate carries its Monte Carlo error. `power` and `power_mixed` gain `power_mcse`,
  `power_lo` and `power_hi`, the standard error and 95% Wilson interval of power. Both also gain
  `type_m_mcse`, the standard error of Type M, and a Wilson interval for Type S, `type_s_lo` and
  `type_s_hi`, since the standard error of Type S is zero whenever no sign is wrong. Each
  `power_curve` record gains the same keys and `type_s`, which it lacked. Type S and Type M rest
  on the significant replicates alone, which are often a handful, yet they were reported without
  any uncertainty. The values are computed as in the R twin, operation for operation, so the same
  counts give the same doubles in both languages, and both test suites pin them. The power guide
  draws the band around its curve from `power_lo` and `power_hi`, where it had computed a Wald
  band by hand, which narrows to nothing where power reaches 0 or 1.
- `power`, `power_curve` and `power_mixed` gain `alternative`, which plans the one-sided test of
  a directional prediction with `"greater"` or `"less"`. Every earlier result was two-sided, and
  the default, `"two.sided"`, leaves every number as it was. The values are spelt as in the R
  twin, and scipy's `"two-sided"` is refused, so that a call means the same in both languages.
  The one-sided p-value comes from the statistic behind the two-sided one, the t of
  `scipy.stats.ttest_ind` or the Wald z of MixedLM. Type S and Type M are taken over the
  replicates significant in the tested direction. `power` records the test in a new key,
  `alternative`, and the `backend` of `power_mixed` names the direction. On the same replicates,
  a one-sided test at alpha rejects exactly where the two-sided test at 2 alpha rejects with the
  estimate in the tested direction. A two-sided result at 2 alpha from an earlier version
  therefore gives the one-sided power at alpha, in the direction of the true effect, as
  `power * (1 - type_s)`.
- `power_mixed` fits the model the specification implies, as the R twin does. `model_formula`
  and `model_data`, ports of the R functions, write that model in lme4's notation and build the
  columns it reads. Every coefficient is a fixed effect, an interaction `"a:b"` entering as its
  product column `a_b`, and every random intercept and slope the specification declares is a
  variance component, those of an extra grouping factor such as a site included. Up to 0.3.0 the
  model was `yv ~ cc` whatever the specification declared, with by-subject and by-item intercept
  and slope components and a test of the first contrast of a single within factor, and any other
  design was refused. A 2 × 2 within design raised `NotImplementedError`, and on
  `reading_time_continuous` six of the seven coefficients, the three predictors and the two
  by-subject slopes were left out. Both now run, and so does any design that declares a random
  effect, with or without an item unit. `tools/parity/formula_cross.py` gives the ports the data
  R simulated for every shipped example and parity case and requires R's formulas character for
  character and R's contrast and product columns bit for bit. The parity job in CI runs it.
  A replicate of `reading_time_continuous` takes about 7 s with NumPy held to one thread. The
  random-effects design is held in sparse matrices from 100,000 entries, rows times random
  coefficients, which made a replicate of the three largest crossed examples 1.4 to 4.5 times
  quicker, and in dense ones below, where sparse ones were about five times slower.
- `power_mixed` gains `focal`, which chooses the effects to test, as in the R twin. By default
  it takes every coefficient, with the true value the specification gives it. A list of names
  leaves the true values unknown, and a dict from names to numbers gives them. Each effect's
  results sit under a new key, `effects`, one entry per focal effect with the twelve keys from
  `power` to `type_m_mcse` that the result held for its one effect. A focal name that the model
  lacks draws the R twin's warning, word for word.
- `power_report` writes a sample-size justification in Markdown from a power result and,
  optionally, a solved curve. The guides called the sample size from `target_n` the number a
  preregistration quotes, but the justification behind it had to be assembled by hand. The report
  gives the items Lakens (2022, <https://doi.org/10.1525/collabra.33267>, Table 4) asks for. It
  states the primary analysis with its test, alpha and sidedness, and the desired power with the
  solved sample size and its interval. Each focal effect follows with the justification given for
  its size, and `null_result` gives the power when every focal effect is zero. The simulation is
  accounted for as Morris, White and Crowther (2019, <https://doi.org/10.1002/sim.8086>)
  recommend, with the replicates, the fit counts and the Monte Carlo standard error and Wilson
  interval of every rate, computed in the report from the counts. The report ends with the
  software versions, the specification's seed and an MD5 checksum of the specification file's
  bytes. A `power_mixed` result is read from `effects`, so
  its deprecated keys are not touched. The R twin's `power_report()` writes the same text, byte
  for byte, from the same counts, which `tools/parity/report_cross.py` checks in CI.

### Changed

- The `power_mixed` docstring and the power guide say that a preregistered simplification rule,
  which moves to a simpler model whenever a fit is singular, is simulated only in the R package.
  There, `power_mixed()` and `precision_design()` take a ladder of models and report the share of
  the replicates analysed with each. `power_mixed` here analyses every replicate with its one
  model.
- `power_curve` gives each grid point its own replicates. Every point reused the specification's
  seed, so neighbouring points shared their random numbers, while `solve_curve` treats them as
  independent. Over 1,000 emulated curves a setting, its nominal 95% intervals covered the true
  sample size in 88% of solves on a tight grid. Point j of a sweep over J sizes now takes
  `replicate_seeds(seed, J)[j]`, the rule of the R twin's `sweep_spec()`, which brings that
  figure to 96.5%, so every number `power_curve` produces changes. `seeds="shared"` keeps the
  old behaviour. The coverage figures quoted for `solve_curve` come from a new calibration over
  1,000 curves a setting, where the earlier 35 in 36 rested on 12 curves whose points shared
  their seeds. Rerun with independent points, those 12 curves no longer favour the probit link
  over the logit, and the documentation now says so.
- The home page gives the Zenodo concept DOI as a full link,
  <https://doi.org/10.5281/zenodo.21266313>, in place of the Zenodo badge, and the About
  page links every DOI the same way, so its reference line can be clicked. The BibTeX
  entry's `url`, on the page and in `docs/pilotr.bib`, gives the PyPI page, so a style that
  prints both fields no longer shows the DOI twice.
- `power` refuses designs whose rows are correlated, with the R twin's message. It t-tested
  every row as an independent observation and accepted any Gaussian design with a two-level
  between factor, including designs crossed with items, designs with a within factor and
  subjects nested in sites. Take 30 subjects crossed with 20 items, with by-subject and
  residual standard deviations of 1 and a by-item one of 0.3. With no true effect, the row-wise
  t-test was significant in 105 of 200 replicates, against 16 for a t-test of the subject means
  on the same data. It now takes one row per subject and no clustering, and otherwise raises
  `NotImplementedError` naming the item unit, within factor or grouping factor it found.
  `power_curve` refuses the same designs.
- The replicates of `power`, `power_curve` and `power_mixed` skip validation, since the
  specification is validated once before them. Each replicate used to validate it again.
- `validate_spec` refuses a specification with two or more factors between the same unit unless
  it declares specification version 0.4, which allocates such factors jointly. pilotr assigned
  the levels of each between factor to blocks of units on its own, so two such factors fell into
  the same or overlapping blocks. A 2 × 2 between-subjects design over 40 subjects produced cells
  of 20, 0, 0 and 20, and the second effect and the interaction could not be estimated. The
  refusal names the encoding that works in every version, one between factor whose levels are
  the cells, and reads the same in the R twin.
- The crossed mixed-effects example on the power page states that every subject sees every item
  in both conditions. It points a counterbalanced study, which has half the observations and
  lower power, to the `counterbalance` field of specification version 0.4, which `power_mixed`
  analyses. For earlier releases it points to the two-list encoding of the specification, which
  `power_mixed` analyses as well.
- The specification page states how units are allocated to levels and clusters, and gives
  worked encodings for a 2 × 2 between design, a two-list counterbalanced design and
  randomisation within clusters. It also says that simulated data are complete, so expected
  attrition is allowed for by recruiting N / (1 − p).
- `power_mixed` fits only the random terms the specification declares, each as an independent
  variance component, and warns once per call when the specification correlates a group's random
  effects, naming the groups whose correlations the model leaves out. A design with random
  intercepts alone is no longer fitted with slope components as well, so its results change. On
  the 12 × 8 design of the power guide, which declares by-subject and by-item intercepts and
  slopes, every number is as before. Its documentation states what the model fits and what it
  leaves out, where the project README had put the two languages at capability parity.
- `power_mixed` refuses a design without random effects with `NotImplementedError`, which
  names `power` for the two-group design and the R twin's `lm()` and `glm()` fits for any
  other. It used to refuse every design without an item unit. A contrast without a coefficient is
  not in the model, as in the R twin, so a null effect is planned with a coefficient of 0, and a
  specification without coefficients is refused in the R twin's words. Up to 0.3.0, `{}` tested
  the contrast against a true effect of 0.
- `n_returned` from `power_mixed` counts every replicate that returned a fit, and each effect's
  power is taken over the fits that give that effect a finite estimate and p-value, as in the R
  twin. A fit that left the tested contrast without a p-value used to count as no fit.
- For a `bernoulli`, `poisson`, `ordinal` or `beta` response, `power_mixed` withholds the
  Type S of an interaction, as the R twin does. A bounded response can reverse the sign of a
  difference of differences between the response scale, where the linear model estimates it, and
  the link scale of its true value.
- `n_converged` from `power_mixed` counts the fits that converged without a warning and without
  a singular variance component, as in the R twin. It counted every fit that returned, which is
  now `n_returned`. `n_attempted`, `n_singular` and `n_warning` are new, and power is still taken
  over the returned fits. A fit is singular when a variance component's standard deviation is
  below 1e-4 of the residual one, lme4's test. statsmodels' notice that a fit may be on the
  boundary counts against no fit, since every fit on the log scale draws it. Its notice of a
  Hessian that is not positive definite counts only against a fit that is not singular, since a
  component at zero can draw it by itself. On the 12 × 8 design of the power guide, the five
  counts equal those of the R twin's `power_mixed()` given the same model.
- The `power_mixed` backend names its test, a Wald z. With few subjects or items, a Wald z
  rejects more readily than the R twin's Satterthwaite tests, so in small designs Python's power
  tends to exceed R's on the same data: 6 of 12 replicates against 4 of 12 on the 12 × 8 design
  of the power guide. With 8 subjects, 6 items, random intercepts alone and no true effect, it
  rejected in 0.070 of 300 replicates (Monte Carlo standard error 0.015), where the R twin's
  Satterthwaite test of the same model rejected in 0.057 on the same data.
  `tools/calibration/wald_null_rate.py` measures both.
- `power_mixed` fits `bernoulli`, `poisson`, `ordinal` and `beta` responses with a linear model
  on the response scale. It now warns for them and withholds the mean estimate and Type M,
  returning both as `nan`, since they compared that scale with a true value on the link scale. A
  linear model of a 0/1 response estimates a difference in probability, which says nothing about
  how far a logit-scale effect is exaggerated. Type S, which compares signs alone, and power are
  kept. The warning reads as the R twin's, whose `power_mixed()` now fits `bernoulli` and
  `poisson` responses with `lme4::glmer()` on their link scale.
- The power guide sets each response family against the model each `power_mixed` fits, and
  against the R package's `precision_design()`.
- NumPy integers and floats, `pathlib` paths and pandas data frames are accepted wherever a
  number, a path or a curve is, and `validate_spec` accepts a path, as in R. A NumPy integer was
  refused as not a whole number, so `power_curve` refused a grid made with `np.arange`, and a
  path object was refused as not a JSON object. `workers` takes any whole number, as in R, and
  `solve_curve` reads a NumPy value of `n` as a replicate count, where it looked for a column of
  that name. `power_curve` records `n_subject` as a Python int, so a curve swept over a NumPy
  grid can be written as JSON.
- The cross-language guarantee is stated as it holds, in the README, the documentation site and
  the specification. Python and R produce bit-identical data wherever they share a maths library,
  as on Linux. On Windows they do not, and a value that passes through `exp()`, `log()` or
  `pow()` can differ in its last bits. Measured there, about 1.4e-5 of normal draws differed,
  since the tail branch of AS 241 calls `log()`, so Gaussian data are not exempt. Lognormal and
  ex-Gaussian responses differed by one unit in the last place and unrounded Beta responses by up
  to 35, and a discrete outcome can in principle move by a whole category. The guarantee used to
  cover the Gaussian family everywhere and every family with `round` set, although `round` is
  refused for the discrete families. Rounding has removed every continuous difference observed,
  but it guarantees nothing. One Beta case can differ on Linux too, under R 4.6.0 or later. From
  that version, R's `^` computes a cube or a fourth power by multiplication, where Python calls
  `pow()`, so a Beta response whose Gamma shape is exactly 1/3 or 1/4 can differ in its last bits.
  Within one language, the data are identical on one platform, and across operating systems on
  the same terms: on Windows, CPython's `exp()` departed from the correctly rounded value for
  0.48% of 700,000 arguments.
- The specification documents every sampling rule a third implementation needs, in a section of
  pseudocode: the partial Fisher-Yates shuffle, the block rules, uniform predictors and
  reliability, the Cholesky factorisation, the two-branch inverse logit, each family's response,
  the Poisson walk, the ordinal rule, the Marsaglia-Tsang Gamma sampler and the Beta ratio. It
  corrects the Gamma sampler's draw consumption. A proposal draws a uniform only when its normal
  gives a positive `v`, and a boosted variate draws its uniform after the Gamma(shape + 1)
  variate. The schema refuses an unknown top-level field and a response parameter its family
  does not use, requires each family's parameters, and accepts the single-value forms of
  `vary_within` and `thresholds`, as `validate_spec` does in strict mode. The parity harness
  checks that the schema accepts every specification of its battery that both twins accept.
- The `replicate_seeds` docstring gives the right reason why consecutive seeds correlate. Seeding
  is affine in the seed, so the first draws of consecutive seeds step by a constant mod 1 and
  correlate 0.95, and a longer warm-up changes only the step (Matsumoto et al., 2007,
  <https://doi.org/10.1145/1276927.1276928>). It said that consecutive replicates began a few
  steps apart in one sequence.
- The `solve_curve` docstring says that when an effect varies from item to item, power at a
  fixed number of items levels off below 1, however many subjects are added, as Westfall, Kenny
  and Judd (2014, <https://doi.org/10.1037/xge0000014>) and Judd, Westfall and Kenny (2017,
  <https://doi.org/10.1146/annurev-psych-122414-033702>) show. A curve over subjects that levels
  off near the target can then be solved far from the true crossing, or as though it crossed when
  it never does. The docstring points to curves over the number of items, one at each of a few
  subject counts large enough to keep their own ceiling well above the target. The R twin's
  `?solve_curve` says the same.
- The contributing guide names the checks in `tools/parity/` that CI runs on the cross-language
  contract. It named `examples/parity_check.py`, which compares CSV files to a tolerance of 1e-6.
  The guides of both no-code apps say that `power_mixed` fits with statsmodels. They said lme4.
  The theme toggle on the website's landing and not-found pages, which lead to this package's
  documentation, now works when the browser refuses the page any storage.
- The docstrings of `load_spec` and `simulate` say that a `validate` other than `True` or
  `False`, `None` and `0` among them, validates leniently, as in the R twin. `load_spec` named
  truthy values only, and `simulate` said nothing of it.
- The worked examples, the home page, Get started and the `pilotr_example` docstring say that eight
  example specifications ship with the package, and the worked examples and the response-families
  guide that they cover six of the eight families, the guide naming the two left out, `bernoulli`
  and `exgaussian`. They said one per family. The worked examples describe `reading_time_continuous`
  as having continuous item and subject predictors and a within factor, where they named one
  continuous predictor, and the `pilotr_example` docstring names the `ValueError` an unknown name
  raises.
- The README and Get started say that the `power` extra serves `power` and `power_curve`, and the
  `mixed` extra `power_mixed` and `model_data`, where they called the first a power demo and the
  second crossed mixed-effects power. The home page says that the extras serve the power tools and
  `model_data`, where it named the power tools alone. The README, the home page and Get started
  name the rare Beta case in which the twins' data can differ on Linux too, which the
  cross-language guide describes.
- The README and the home page give one list of what the R package alone offers: precision and ROPE
  design analysis, power curves for mixed-effects designs, calibration and sweeps, the `lme4`
  backend, the `brms` bridge, the analysis-script emitters, the helpers that build specifications
  and the app. Each left out part of it, and the home page said that pilotr estimates precision and
  emits analysis code, which only the R package does. The API reference names precision and
  equivalence among the groups of the R package's reference that have no Python counterpart, and
  says that `power_mixed` fits the model `model_formula` writes with each random effect as an
  independent variance component.
- The docstrings of `power`, `power_curve`, `power_mixed`, `model_formula` and `model_data` say
  that they validate the specification strictly, even one that `load_spec` read leniently. The
  three power functions list an invalid specification and an invalid `workers` among the
  errors they raise, and `power_curve` lists the `NotImplementedError` it raises for a design
  that `power` refuses. The `simulate` docstring lists the errors it raises, among them a
  random-effect covariance that is not positive definite, which validation does not check, and
  says that a seed and its negation give the same data, as the `RNG` docstring now does too.
- The `inv_logit` docstring gives its range as [0, 1], since it rounds to exactly 1 from an
  argument of about 36.7 and underflows to 0 below about -745.1. It said (0, 1).
- The `power_mixed` docstring and the power guide say that the link keeps the sign of a main
  effect unless an interaction reverses that effect between the levels of another factor, as the
  R twin's `?power_mixed` says. They said that the link keeps it.
- The `solve_curve` docstring says that `n` may be a list or tuple of counts, one per row that
  `effect` keeps, recycled if shorter and cut if longer. `target_n` calls a power of 0.8 the
  conventional target, where it named the plots of this package, which draws none.
- The power guide says that each entry of `res["effects"]` holds the twelve keys from `power` to
  `type_m_mcse`, where it said every key of a `power` result, and gives the formula of the model
  `power_mixed` fits, `model_formula(spec_mixed, structure="zcp")`, for the R package's
  `power_mixed()` to reproduce its fit counts.
- The cross-language guide says that the parity harness allows the unrounded dumps of
  `crossed_mixed_rt`, `reading_time_continuous` and `stress_exgaussian` up to 8 units in the last
  place on every platform, Linux included, and holds every other case to the last bit. It said
  that every case was held to the last bit. It says that `examples/parity_check.py` compares the
  CSV files both demo scripts write for six of the eight example designs, all but the ordinal and
  Poisson ones, and that its difference below 1e-14 lies in the unrounded predictors of
  `reading_time_continuous`, whose rounded response agrees exactly. It said that the check covered
  every worked example and that the continuous design was unrounded. It says that
  `tools/parity/solve_cross.py` compares the text of a selection of the refusals of `solve_curve`,
  where it said every refusal. It also says that the crossed reaction-time example has one row per
  subject, item and condition, where it said one per subject and item.
- The response-families guide says that an ex-Gaussian specification has to declare a
  `spec_version` of 0.3 or later. It named 0.3 alone.
- The `RNG` and `as241` docstrings give the references they cite, L'Ecuyer (1988,
  <https://doi.org/10.1145/62959.62969>) and Wichura (1988, <https://doi.org/10.2307/2347330>),
  and the `power` docstring and the power guide give the DOI of Gelman and Carlin (2014,
  <https://doi.org/10.1177/1745691614551642>). The `solve_curve` docstring cites Fieller (1954,
  <https://doi.org/10.1111/j.2517-6161.1954.tb00159.x>), whom its references listed without a
  citation: in the checks against R's `stats::power.t.test`, Fieller's interval covered the
  analytic answer as often as the delta-method interval, at nearly the same width.
- The specification says that base R's `sum()`, which rounds its long-double total to double a
  second time, can differ from a plain double fold for an inner product of two terms, where it said
  three terms or more. It says that `build_spec()` rounds an ex-Gaussian response to four decimals,
  as it rounds the Gaussian and lognormal ones, and it sets out how a sweep seeds its points under
  Replicate seeds, beside the rule that paragraph extends. Where it described what earlier releases
  did, it gives the rule and its reason, and it no longer lists Matuschek et al. (2017,
  <https://doi.org/10.1016/j.jml.2017.01.001>), whom it does not cite. The schema's comment on the
  family blocks names the tables of both validators, where it gave the Python names for both.
- The repository's README says that `workers` spreads the replicates over local worker
  processes, where it said in-process, and that both packages ship the example specifications
  `pilotr_example()` lists. It names the packages the R package suggests for its app, its
  mixed-effects power and its precision analysis. The note on advanced designs in the repository
  says that `power_mixed` has no joint test of several coefficients and refuses a design without
  random effects, where it said that neither was planned, and passes `effect` to the R package's
  `solve_curve()` for a curve over several focal effects.
- The contributing guide lists all eight checks of the parity job and says how to build this
  package's documentation site, where the two builds of the no-code app live and how to run
  their tests. The website's landing page says that the browser app estimates power for a
  two-group Gaussian design with one row per subject, where it called the app the same workflow
  as the packages, and the security policy says that only that build runs client-side.
- The repository's README and its note on advanced designs say that the parity harness allows
  the unrounded dumps of `crossed_mixed_rt`, `reading_time_continuous` and `stress_exgaussian` to
  differ by up to 8 units in the last place, on Linux as on any other platform, and holds every
  other case to the last bit. They said that continuous integration held every case to the last
  bit. The README's note on R and Python coverage, the contributing guide and the website's
  landing page name the rare Beta case in which the twins' data can differ on Linux too.
- The power guide says that this package builds no sweep over items or effect sizes and no
  curve for a mixed-effects design, that the R package's `power_curve_mixed()` and `sweep_spec()`
  build them, and that `solve_curve` reads such a curve once R has written it to a CSV file and
  pandas has read it back.

### Fixed

- The citation on the About page names the Python package version alone. It read "R and
  Python package version", although the two packages are released separately and their
  version numbers can differ. The BibTeX file the page offers for download says the same.
- The 0.3.0 entry carries the date the release reached PyPI, 2026-08-27.
- `load_spec` reads a specification as R's jsonlite does. A one-element array where a single
  value belongs, `[]` for an empty `coefficients`, `slopes` or `correlations` object (both as
  `jsonlite::write_json()` writes them by default) and a JSON `null` in an optional field are
  read as the value, the empty object and an absent field. Each previously raised or was
  refused.
- A whole number written as a float, such as `"n": 64.0`, no longer raises `TypeError` in
  `simulate`. The schema counts it as an integer and R simulates it.
- `validate_spec` raises `ValueError` for every malformed shape, as documented. A list-valued
  family raised `TypeError` and a string-valued `response` raised `AttributeError`.
- `validate_spec` refuses `"item": null` among the units and a `random` of `0` or `false`, in
  R's words. The first passed validation and then raised `TypeError` in `simulate`, and the
  second was read as no random effects.
- `load_spec` reads files as UTF-8 and accepts a byte-order mark, as R does. It used the locale
  encoding, so on Windows the level "fácil" loaded as "fÃ¡cil" and "Łatwy" raised
  `UnicodeDecodeError`.
- `load_spec` refuses a repeated key within one object, and `validate_spec` a seed beyond
  ±(2^53 − 1). Both used to give data that differed from R's.
- `validate_spec` refuses specifications whose names collide. Each of these validated and then
  moved an effect, rescaled a variance or overwrote a column without a word: two columns with one
  name (a factor called `subject`, a response named like a factor), a contrast column defined by
  two factors or named like a predictor, a unit, a grouping factor, another factor or the
  response, an interaction whose analysis column `a_b` already exists, a level listed twice, a
  correlation pairing a term with itself or giving one pair twice, a factor both between and
  within a unit, and a `random.item` entry in a design without items, which was dropped. Where a
  column was overwritten, `simulate` wrote a second column under the same name while R replaced
  the first, so the twins exported different tables. Each refusal reads the same in the R twin.
- `validate_spec` refuses a factor, predictor or extra grouping-factor name that is empty or made
  only of spaces, and a response name made only of spaces (an empty one was already refused).
  Such a name validated and then became a column with no visible name.
- `power` reads a missing coefficient for the tested contrast as a true effect of 0, with Type S
  and Type M undefined. It raised `KeyError`, as `power_mixed` did, which now tests the
  coefficients the specification lists (see Added).
- `power_mixed` analyses a `lognormal` response on the log scale, as it already did for
  `shifted_lognormal` and as the R twin does. It fitted the raw values, so the mean estimate came
  out in the response's units. Take 12 subjects crossed with 8 items and a true effect of 0.1 on
  the log scale. Six replicates gave a mean estimate of 53.1 and a Type M of 553, where the log
  scale gives 0.136 and 1.36, as R does.
- `power_mixed` returns no estimate for a replicate whose response has no logarithm, as the R
  twin does, and `n_returned` leaves that replicate out. Such a response is a
  `shifted_lognormal` value rounded down to the shift, or now a `lognormal` value rounded to 0.
  The shifted family raised `ValueError: math domain error`, which ended the whole call.
- `power_mixed` fits with Powell's method followed by L-BFGS. statsmodels' default chain of
  optimisers stopped short of the REML optimum in most crossed fits, which the documentation had
  described as statsmodels overstating random-slope variance. Over the first 12 replicates of the
  crossed reaction-time example, the standard error of the effect averaged 0.0278 against lme4's
  0.0190, and 2 replicates were significant against lme4's 6. The standard errors now match those
  of lme4's uncorrelated model to four decimal places in 11 of the 12, and the same 6 replicates
  are significant. A replicate whose Powell fit raises is fitted again with the default chain
  and counts as a fit with a warning.
- A rounded response is rounded by R's algorithm, so the twins agree at decimal ties. CPython's
  `round()` rounds the exact binary value, while R takes whichever of the two decimal candidates
  is nearer in double precision, so `round(1.05, 1)` gave 1.1 in Python and 1 in R. The two
  disagreed on 166 of 5,000 values probed within a few units in the last place of a tie and on
  none of 1,000 random values. The shipped examples simulate the same data as before. With
  `validate=False`, a `round` that is not a whole number, which validation refuses, counts as
  `floor(round + 0.5)` decimal places, as in R. It raised `TypeError`.
- A Poisson mean above about 708 is refused. The sampler's starting value `exp(-mean)` is
  subnormal there, and from a mean of about 738 the counts were biased or intermittently refused,
  while the specification promised service up to 746. Over 40,000 draws, a mean of 744 drew
  counts with a mean of 733.5 and a variance of 407, and a mean of 742 was refused for 4% of
  them. The refusal reads as the R twin's and asks for a mean of at most 708.
- The Poisson sampler evaluates each term of its recurrence as `(p * lam) / k`, as R does. The
  other bracketing, `p * (lam / k)`, lands on a different double for most cumulative sums, and a
  uniform falling between the two gave a different count. At a mean of 3.7, a uniform of
  0.4941532441504185 drew 4 here and 3 in R. In 300,000 draws from the generator, at means from
  0.5 to 100, no count differed.
- `Dataset.to_csv` writes UTF-8 and quotes a field that holds a comma, a double quote or a
  newline. It wrote the locale encoding, so on Windows a Polish label raised
  `UnicodeEncodeError` and a Spanish one came out in bytes that no UTF-8 reader accepts. A comma
  in a label or a column name shifted every later field of its row. Lines still end in LF on
  every platform, and a data set without such characters is written byte for byte as before.
- A NumPy float64 in a specification no longer reaches the simulated data as a NumPy scalar,
  which `to_csv` wrote as `np.float64(...)` under NumPy 2. Every number in the specification
  becomes a Python int or float before anything is drawn, on a copy of the caller's.
- A specification written by the R twin's `spec_json()` carries exactly the numbers R held, so
  `load_spec` reads the same doubles. The R writer checked its shorter forms with a reader that
  rounds twice on R builds with long doubles, and about one number in 5,000 to 7,000 was written
  one unit in the last place away, which `json` then read as written. The Python package is
  unchanged. A fixture written by R and a parity check over 20,000 numbers now hold the two
  packages to the same reading.
- Beta responses stay inside (0, 1). Near a ceiling, `x / (x + y)` rounds to exactly 1 in double
  precision (161 of 4,000 rows at an intercept of 3 and `phi` 2), and rounding sent low
  proportions to exactly 0 (17% of the rows at three decimals in one realistic design). brms's
  `Beta()` refuses either value. An exact 1 or 0 now becomes the nearest double inside the
  interval, a rounded value at 0 or 1 becomes the nearest value inside it at that many decimals,
  and `validate_spec` refuses `round` 0 for Beta, as the R twin does. No other value changes, so
  the shipped examples simulate the same data as before. The sampler is
  `pilotr.core.beta_response(rng, eta, phi)`, which replaces `beta_draw(rng, a, b)`, a function
  the documented API never included.
- Extreme linear predictors are refused with the R twin's message. `simulate` raised a bare
  `OverflowError` for a lognormal or shifted lognormal response whose `exp(eta + sigma * z)`
  overflowed, and `ZeroDivisionError` for a Beta response whose shape parameter was zero or whose
  two Gamma variates both underflowed. R returned `Inf`, 0, 1 or `NaN` in the same cases.
- The install commands in the README and the getting-started guide carry no trailing comments,
  so they paste into any shell. At an interactive prompt, zsh, the default shell on macOS, reads
  `#` as the start of a comment only when told to, so `pip install pilotr  # core engine` handed
  the `#` and the words after it to pip, which stopped with "Invalid requirement: '#'". The
  blocks now hold the bare commands, and the prose before them says what each extra adds. The
  extras were already in double quotes, which keep zsh from reading the brackets as a filename
  pattern. A test fails on a `pip install` line in a shell block that carries a `#` comment, and
  on an extra outside double quotes, in the README, the guides, the examples and the sources. In
  a checkout, it also reads the repository's README and both copies of the contributing guide.
- `power` refuses a design with fewer than two subjects at a level of the between factor, so
  fewer than 4 in all, in the R twin's words, and `power_curve` refuses such a grid point. With
  one subject at each level, `power` returned a power of 0, since scipy answers a group without
  variance with a nan p-value and a warning, while the R twin stopped in `t.test()`. With one
  subject in all, it raised `StatisticsError`. The error is a `ValueError`.
- A Beta response whose Gamma shape is exactly 1/2 is drawn as the R twin draws it, which changes
  this package's data for such designs. The Gamma sampler boosts a shape below 1 by
  `u ** (1 / shape)`, and Python's `**` calls `pow()`, where R's `^` computes a square as `u * u`.
  On Linux, the two differed for 0.08% of the generator's uniforms. A Beta design with a linear
  predictor of 0 and a `phi` of 1, whose shapes are both 1/2, therefore differed from R in 13 to
  19 of 20,000 rows at each of four seeds, by up to 2 units in the last place. The boost is now
  `u * u` when `1 / shape` is exactly 2, and those rows hold R's values. No other row changes, and
  a rounded response changes only where the difference crossed a rounding boundary. The boost for
  a shape of exactly 1/3 or 1/4 still calls `pow()`, as R did before 4.6.0, and can differ from
  R 4.6.0 or later.
- `validate_spec` applies the checks on the unit columns to an empty `units` object, as the R twin
  does, so the two report the same problems for one file. R reads `{}` as a named empty list and
  checks it as an object, while this package skipped those checks for an empty object. A
  specification with `"units": {}`, a within factor named `subject` and a `random.item` entry was
  refused here for the missing `units.subject` alone. R also reported that the name `subject` is
  used for more than one column and that `random.item` describes an item unit the design does not
  have, and both twins now report all three. An empty array for `units.subject` or `units.item`,
  which R reads as an empty list, is read as an empty object too, so both twins report its
  missing `n`, where this package said that the unit must be an object.

### Deprecated

- A within factor whose `vary_within` omits a unit of the design now warns, in strict and lenient
  validation alike. pilotr has always crossed a within factor with every unit, so the list
  changed nothing. From spec version 0.4, such a list is refused.
- The top-level keys of the result of `power_mixed` that held its one effect's quantities,
  `power` to `type_m_mcse`, are kept for this release beside `effects` and will be removed in
  the next. They hold the entry of the first contrast of a single within factor, as before, and
  reading one raises a `FutureWarning` that names its place in `effects`.

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

[Unreleased]: https://github.com/pablobernabeu/pilotr/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/pablobernabeu/pilotr/releases/tag/v0.3.0
[0.2.1]: https://github.com/pablobernabeu/pilotr/releases/tag/v0.2.1
[0.2.0]: https://github.com/pablobernabeu/pilotr/releases/tag/v0.2.0
[0.1.0]: https://github.com/pablobernabeu/pilotr/releases/tag/v0.1.0
