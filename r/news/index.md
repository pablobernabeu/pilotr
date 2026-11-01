# Changelog

## pilotr (development version)

### Upgrading from the previous release (0.3.1)

#### Changed results

- [`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md),
  [`precision_design()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md)
  and their curves fit `bernoulli` and `poisson` responses with
  [`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html), where a
  linear model compared response-scale estimates with link-scale true
  values. Estimates, interval widths, Type M, ROPE decisions and power
  change. Where `formula` is taken, `formula = model_formula(spec)`
  restores the old model.

- Each point of a sweep draws its own replicates, so the results at
  every point of
  [`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md),
  [`power_curve_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_curve_mixed.md),
  [`precision_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_curve.md)
  and the no-code app’s power curve change. Shared seeds made
  [`solve_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/solve_curve.md)
  intervals too narrow. `seeds = "shared"` (`.seeds` in
  [`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md))
  restores the 0.3.1 seeding.

- For `ordinal` and `beta` responses,
  [`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md),
  [`precision_design()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md)
  and their curves warn and withhold the mean estimate, Type M, ROPE
  decisions and an interaction’s Type S, where they compared a
  response-scale estimate with a link-scale true value or ROPE.
  `formula = model_formula(spec)`, where taken, restores them.

- [`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md)
  merges each condition from
  [`design_conditions()`](https://pablobernabeu.github.io/pilotr/r/reference/design_conditions.md)
  into the specification’s coefficients, where it replaced them and
  dropped every coefficient a condition did not name.
  [`unclass()`](https://rdrr.io/r/base/class.html) on the conditions
  restores the replacement.

- For a counterbalanced pilot,
  [`spec_from_model()`](https://pablobernabeu.github.io/pilotr/r/reference/spec_from_model.md)
  sets `counterbalance` and declares specification version 0.4, where it
  crossed every subject with every item in every condition and so
  overstated power. `counterbalance = FALSE` restores the crossed
  specification, with a warning.

- For two or more factors between one unit,
  [`spec_from_model()`](https://pablobernabeu.github.io/pilotr/r/reference/spec_from_model.md)
  builds one factor whose levels are the cells, where its specification
  left cells empty, half of them in a 2 × 2 pilot.

- With `validate = FALSE`,
  [`simulate_design()`](https://pablobernabeu.github.io/pilotr/r/reference/simulate_design.md)
  allocates two or more factors between one unit jointly, one block per
  combination of levels, where each factor took its own blocks, which
  coincided or overlapped and left cells empty.

- [`model_data()`](https://pablobernabeu.github.io/pilotr/r/reference/model_data.md)
  fills every contrast column of a factor whose first contrast column
  carries the factor’s own name, where it left the later columns `NA`,
  so a model using them could not be fitted and every replicate of
  [`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md)
  and
  [`precision_design()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md)
  failed.

- `calibrate_response(tune = "all")` holds the intercept, where scaling
  it raised the responses of both lognormal families, less any shift, to
  the power of the scaling factor. For `poisson` and `beta`, the scaling
  factor and the attainable targets change too. `scale_intercept = TRUE`
  scales it as before.

- [`brms_bridge()`](https://pablobernabeu.github.io/pilotr/r/reference/brms_bridge.md)
  and
  [`generate_design_analysis()`](https://pablobernabeu.github.io/pilotr/r/reference/generate_design_analysis.md)
  leave the intercept to brms’s default prior, centred on the response
  for the continuous families and `poisson` and on zero for the others,
  where a fixed `normal(0, 2.5)` prior pulled the intercept of a
  response far from zero towards zero.

- For `gaussian`, `lognormal`, `shifted_lognormal` and `exgaussian`
  responses,
  [`brms_bridge()`](https://pablobernabeu.github.io/pilotr/r/reference/brms_bridge.md)
  and
  [`generate_design_analysis()`](https://pablobernabeu.github.io/pilotr/r/reference/generate_design_analysis.md)
  scale the priors to the response and to each column’s standard
  deviation, where a unit scale ignored the data’s scale. Posteriors and
  Bayes factors change.

- Specification fields are read by their exact names, where partial
  matching let a random entry `subject_site` act as `subject` and, under
  lenient validation, a misspelt `rounding` act as `round`.

- `rope = c(lower, upper)` in
  [`precision_design()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md)
  and
  [`precision_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_curve.md)
  is an interval, where each value was recycled across replicates. A
  single non-negative half-width gives unchanged results. The
  repository’s HPC precision-array job reads a `ROPE` such as
  `-0.02,0.05` as an interval, where it gave `NA`.

- The no-code app’s power curve extends its grid until it brackets the
  0.80 target, where it stopped at twice N and so was often refused or
  solved too low.

- For a design its Power tab cannot analyse, the no-code app recommends
  only analyses that fit it, where it recommended the
  [`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md)
  of both packages for every such design. For ordinal and Beta designs,
  it recommends
  [`generate_design_analysis()`](https://pablobernabeu.github.io/pilotr/r/reference/generate_design_analysis.md)
  and
  [`brms_bridge()`](https://pablobernabeu.github.io/pilotr/r/reference/brms_bridge.md).

- Scripts from
  [`generate_r_script()`](https://pablobernabeu.github.io/pilotr/r/reference/generate_r_script.md)
  recommend
  [`power_design()`](https://pablobernabeu.github.io/pilotr/r/reference/power_design.md)
  for a two-group Gaussian design with one row per subject,
  [`generate_design_analysis()`](https://pablobernabeu.github.io/pilotr/r/reference/generate_design_analysis.md)
  for ordinal and Beta responses and
  [`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md)
  for any other design, where every script recommended
  [`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md).

- The repository’s HPC precision-array job seeds each task as
  [`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md)
  seeds its point, where every task took the specification’s seed, which
  made
  [`solve_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/solve_curve.md)
  intervals too narrow. It counts fits as
  [`precision_design()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md)
  does, where it counted every fit that returned an estimate as
  converged.

- `n_warning` no longer counts boundary-singular fits, which
  `n_singular` counts, so it shows whether any fit failed to converge.

- A Beta response that came out as exactly 0 or 1, which brms’s
  [`Beta()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)
  refuses, becomes the nearest double inside (0, 1) or, for a rounded
  response, the nearest value inside it at its number of decimals. No
  other value changes.

- [`spec_json()`](https://pablobernabeu.github.io/pilotr/r/reference/spec_json.md)
  writes some 15- and 16-digit numbers differently, so that R and Python
  read back exactly the double R holds, where, on R builds with long
  doubles, some came back one unit in the last place away.

- Scripts from
  [`generate_r_script()`](https://pablobernabeu.github.io/pilotr/r/reference/generate_r_script.md)
  and
  [`generate_design_analysis()`](https://pablobernabeu.github.io/pilotr/r/reference/generate_design_analysis.md)
  embed the specification as JSON read with
  [`load_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/load_spec.md),
  where R builds without long doubles, such as on Apple silicon, could
  parse a number in their R code one unit in the last place away and
  simulate different data.

#### Changed interface

- [`power_design()`](https://pablobernabeu.github.io/pilotr/r/reference/power_design.md)
  refuses a design whose rows are correlated by an item unit, a within
  factor or a grouping factor besides `subject`, where it t-tested every
  row as independent and overstated power. The no-code app applies the
  same rule.

- [`validate_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/validate_spec.md)
  refuses two or more factors between the same unit unless the
  specification declares `"spec_version": "0.4"`, where such factors
  fell into overlapping blocks that left cells empty. One between factor
  whose levels are the cells works in every version.

- [`validate_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/validate_spec.md)
  refuses colliding names, such as a factor called `subject` or a level
  listed twice, which moved an effect or overwrote a column, and blank
  names. It also refuses a `random.item` entry in a design without
  items, which was dropped.

- [`validate_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/validate_spec.md)
  refuses a factor that sets both `vary_within` and `between`, which
  gave each unit one row per level, all at one level. It also refuses a
  correlation that pairs a term with itself, which rescaled its
  variance, and two correlations for one pair.

- [`precision_design()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md)
  and
  [`precision_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_curve.md)
  refuse a `rope` other than a non-negative half-width or an increasing
  pair, where they used a negative half-width as given and recycled any
  other vector. So does the repository’s HPC precision-array job for its
  `ROPE`.

- [`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md)
  refuses a
  [`design_conditions()`](https://pablobernabeu.github.io/pilotr/r/reference/design_conditions.md)
  condition that names a coefficient the specification lacks, where it
  used the condition as the whole coefficient list. With
  [`unclass()`](https://rdrr.io/r/base/class.html) on the conditions, it
  uses each as the whole list, as before.

- [`power_design()`](https://pablobernabeu.github.io/pilotr/r/reference/power_design.md)
  refuses a design with fewer than two subjects at a level of the
  between factor, where three subjects in all left the pooled variance
  one degree of freedom.

- Results gain elements and columns.
  [`power_design()`](https://pablobernabeu.github.io/pilotr/r/reference/power_design.md)
  returns `alternative`, and
  [`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md),
  [`power_curve_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_curve_mixed.md),
  [`precision_design()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md)
  and
  [`precision_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_curve.md)
  return `fitter`, which names the fitter and the test. The power
  functions also return the Monte Carlo errors listed under New
  features.

- In a data frame from
  [`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md)
  over
  [`power_design()`](https://pablobernabeu.github.io/pilotr/r/reference/power_design.md),
  `power` and every later column sit further right, as do `type_m` and
  the fit counts over
  [`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md)
  and from
  [`power_curve_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_curve_mixed.md).
  Code that selects them by position gets other columns.

- The specification
  [`spec_from_model()`](https://pablobernabeu.github.io/pilotr/r/reference/spec_from_model.md)
  returns carries an `effects` attribute, which records how its effects
  were chosen.

- [`simulate_design()`](https://pablobernabeu.github.io/pilotr/r/reference/simulate_design.md)
  refuses a Poisson mean above about 708, where counts were biased or
  intermittently refused from a mean of about 738 upwards.

- [`simulate_design()`](https://pablobernabeu.github.io/pilotr/r/reference/simulate_design.md)
  refuses extreme linear predictors with the same message as the Python
  twin, where it returned `Inf` for an overflowing lognormal or shifted
  lognormal response and 0, 1 or `NaN` for a Beta response.

- [`validate_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/validate_spec.md)
  refuses a seed beyond ±(2^53 − 1), which R and Python can read as
  different numbers, and a `round` of 0 for a Beta response, which
  rounded every response to 0 or 1.

- [`load_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/load_spec.md)
  refuses a key repeated within one object, which R and Python read
  differently, even with `validate = FALSE`, and so does every function
  given a path.
  [`validate_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/validate_spec.md)
  and every function that validates a list refuse such a key in a list
  too.

- [`validate_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/validate_spec.md)
  refuses a `spec_version` given as an array, which was read as its
  first element, and a factor level, a `vary_within` entry or a
  `correlated` flag written as `null`, which became `NA` or `false`.

- [`build_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/build_spec.md)
  leaves a seed, number of subjects or number of items that is not a
  whole number for
  [`validate_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/validate_spec.md)
  to refuse, where it truncated it.

- [`brms_bridge()`](https://pablobernabeu.github.io/pilotr/r/reference/brms_bridge.md)
  validates a specification given as a list. It refuses an interaction
  given two keys, such as `cond:age` and `age:cond`, and a `prior_scale`
  or `interaction_scale` that is not a single positive number.

- `default_response_name("exgaussian")` returns `RT`, where it returned
  `outcome`.

- The repository’s HPC precision-array job runs once its `--account` and
  `PROJECT_DIR` are set, from the clone’s `hpc/` directory with pilotr
  in `$PROJECT_DIR/Rlib`. It no longer charges the author’s account,
  runs a runner kept outside the repository or sources the pilotr that
  `TOOLKIT` names.

#### Changed output

- Printing a
  [`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md)
  result names the fitter, shows Type S and Type M with their Monte
  Carlo errors in a second table and applies `digits` to the interval
  bounds, which had three decimals whatever `digits` said.

- The no-code app’s Verify button runs the exported script and certifies
  it only when its data are identical, where it compared response sums
  within a tolerance without running it. Otherwise, it names the first
  column and row that differ.

- Scripts from
  [`generate_r_script()`](https://pablobernabeu.github.io/pilotr/r/reference/generate_r_script.md)
  embed a specification given as a path, where they wrote the path
  itself, and install pilotr from CRAN, where the installation line said
  ‘once available’. Where they recommend
  [`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md),
  the installation line also names lme4 and lmerTest.

- The browser build of the no-code app gives commands that install
  pilotr from CRAN, with lme4 and lmerTest or brms, and the Python
  package from PyPI, where it installed both from GitHub with a pip
  command that failed in zsh.

- The comments in the scripts
  [`generate_design_analysis()`](https://pablobernabeu.github.io/pilotr/r/reference/generate_design_analysis.md)
  writes match what the scripts do, where they named the wrong output
  directory, tied the R-hat limit to 1.01 whatever `gate` set and left
  the partitions out of the settings to edit.

- Each task of the repository’s HPC precision-array job writes the
  columns of a
  [`precision_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_curve.md)
  row plus `pilotr_version` and `lme4_version`, where it wrote `n_sims`
  and no `n_attempted`, `n_returned`, `n_singular`, `n_warning` or
  `fitter`.

- The no-code app’s Data tab, summary, plot and data download show the
  last simulation and its seed, where the download simulated afresh, the
  installed app’s Data tab showed the Seed box and the plot lost its
  grouping after a factor was renamed.

- The no-code app names each download after the specification it holds,
  where the installed app used the Design name box.

- [`spec_from_model()`](https://pablobernabeu.github.io/pilotr/r/reference/spec_from_model.md)
  says in messages which fixed effects it took from the pilot and which
  factors it combined or counterbalanced, and warns when it returns a
  crossed specification for a counterbalanced pilot.

- For `gaussian`, `lognormal`, `shifted_lognormal` and `exgaussian`
  responses,
  [`brms_bridge()`](https://pablobernabeu.github.io/pilotr/r/reference/brms_bridge.md)
  warns about a coefficient whose column does not vary.

- A within factor whose `vary_within` omits a unit of the design warns
  (see Deprecated).

- [`validate_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/validate_spec.md)
  reports a malformed `response`, `predictors`, `random` or
  `spec_version` in its own words, where some stopped with a base R
  error.

- [`model_formula()`](https://pablobernabeu.github.io/pilotr/r/reference/model_formula.md)
  refuses, with the same message as the Python twin, a specification
  with neither a fixed coefficient nor a random effect, where it stopped
  with a parse error.

- [`validate_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/validate_spec.md)
  words its refusals of a `correlated` flag that is not a boolean and of
  a malformed `response.thresholds` as the Python twin does.

- The code that
  [`brms_bridge()`](https://pablobernabeu.github.io/pilotr/r/reference/brms_bridge.md)
  returns says in a comment that
  `model_data(spec, simulate_design(spec))` adds the contrast columns
  that its formula names.

- `citation("pilotr")` cites the R package, published by CRAN, with the
  DOI that CRAN registered,
  <https://doi.org/10.32614/CRAN.package.pilotr>, and its CRAN page,
  where it cited the R and Python packages with the Zenodo concept DOI,
  <https://doi.org/10.5281/zenodo.21266313>.

### New features

- Specification version 0.4 allocates several between factors on one
  unit jointly, one block per combination of levels with the
  first-listed factor varying slowest. pilotr 0.3.0 and 0.3.1 refuse a
  specification that declares it, and earlier releases do not check the
  version.
  [`build_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/build_spec.md)
  and
  [`spec_from_model()`](https://pablobernabeu.github.io/pilotr/r/reference/spec_from_model.md)
  declare version 0.4 only for a specification that uses a 0.4 feature,
  and 0.3 otherwise, so a specification without one stays readable by
  pilotr 0.3.0 and 0.3.1.

- `counterbalance`, a specification-0.4 field on a within factor, gives
  each subject one condition per item, rotating the lists as in a Latin
  square, and adds a `list` column to the simulated data.
  [`build_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/build_spec.md)
  and the no-code app offer it.

- [`calibrate_response()`](https://pablobernabeu.github.io/pilotr/r/reference/calibrate_response.md)
  gains `denominator = "random"`, which calibrates the random components
  plus the residual, so that, at the default `target_var = 1`, the
  coefficient of a contrast coded `c(-0.5, 0.5)` is Cohen’s d (Westfall,
  Kenny and Judd, 2014, <https://doi.org/10.1037/xge0000014>).

- Every rate carries its Monte Carlo error: the standard error and 95%
  Wilson interval of power (`power_mcse`, `power_lo`, `power_hi`), the
  standard error of Type M (`type_m_mcse`) and a Wilson interval for
  Type S (`type_s_lo`, `type_s_hi`).
  [`power_design()`](https://pablobernabeu.github.io/pilotr/r/reference/power_design.md)
  reports them all, and
  [`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md)
  and
  [`power_curve_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_curve_mixed.md)
  add the last three. The no-code app shows each rate with its error.

- [`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md),
  [`precision_design()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md)
  and
  [`precision_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_curve.md)
  simulate a preregistered simplification rule of the kind Barr et
  al. (2013, <https://doi.org/10.1016/j.jml.2012.11.001>) and Matuschek
  et al. (2017, <https://doi.org/10.1016/j.jml.2017.01.001>) discuss,
  given as a ladder of models in `formula` and a criterion in
  `fallback`, and report the share of the replicates analysed at each
  rung.

- [`model_formula()`](https://pablobernabeu.github.io/pilotr/r/reference/model_formula.md)
  gains `structure`, which builds the maximal, zero-correlation or
  random-intercepts-only model, the rungs of a ladder such as
  `c("maximal", "zcp", "intercepts")`.

- [`spec_from_model()`](https://pablobernabeu.github.io/pilotr/r/reference/spec_from_model.md)
  gains `effects = "safeguard"`, which takes the bound nearer zero of
  each effect’s interval at `safeguard_level`, 60% by default (Perugini,
  Gallucci and Costantini, 2014,
  <https://doi.org/10.1177/1745691614528519>), and `sesoi`, which sets
  named effects to a smallest effect size of interest (Albers and
  Lakens, 2018, <https://doi.org/10.1016/j.jesp.2017.09.004>).

- [`power_design()`](https://pablobernabeu.github.io/pilotr/r/reference/power_design.md),
  [`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md)
  and
  [`power_curve_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_curve_mixed.md)
  gain `alternative`, which plans a one-sided test with `"greater"` or
  `"less"`. The default, `"two.sided"`, is the test of earlier releases.

- `focal` in
  [`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md)
  and
  [`power_curve_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_curve_mixed.md)
  accepts named sets of coefficients, each tested jointly, so the power
  of a multi-degree-of-freedom effect, such as a 3 × 2 interaction, can
  be planned.

- [`power_report()`](https://pablobernabeu.github.io/pilotr/r/reference/power_report.md)
  writes a sample-size justification in Markdown from a power or
  precision result and, optionally, a solved curve, with the items of
  Lakens (2022, <https://doi.org/10.1525/collabra.33267>) and the
  simulation accounting of Morris, White and Crowther (2019,
  <https://doi.org/10.1002/sim.8086>). For a precision result, each
  focal effect has its mean interval width and the probabilities of the
  two decisions about the region of practical equivalence. For a power
  result, the Python twin writes the same text from the same counts.

### Minor improvements

- [`power_curve_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_curve_mixed.md)
  and
  [`precision_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_curve.md)
  start one worker pool for the whole sweep, where each grid point
  started its own.
  [`power_design()`](https://pablobernabeu.github.io/pilotr/r/reference/power_design.md),
  [`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md)
  and
  [`precision_design()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md)
  gain `cl`, which runs them on a cluster the caller has made, and,
  through
  [`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md),
  one such cluster runs every point of a sweep.

- With `async`, the no-code app runs its power curve in a background
  process as well as its power analysis, and, over a sequential future
  plan,
  [`run_app()`](https://pablobernabeu.github.io/pilotr/r/reference/run_app.md)
  starts two workers, where it started one for every core.

- For a design without random effects or continuous predictors, the
  no-code app’s power curve also reports the sample size that
  [`power.t.test()`](https://rdrr.io/r/stats/power.t.test.html) gives
  for a power of 0.80.

- In a browser tab, the browser build of the no-code app is titled
  ‘pilotr: no-code app’ and shows the pilotr icon.

- The no-code app shows a note when its controls have changed since the
  last simulation.

### Fixes

- [`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md)
  and
  [`precision_design()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md)
  fit a design without random effects with
  [`lm()`](https://rdrr.io/r/stats/lm.html) or
  [`glm()`](https://rdrr.io/r/stats/glm.html), where every replicate
  failed, so every between-subjects example can be analysed.

- [`power_design()`](https://pablobernabeu.github.io/pilotr/r/reference/power_design.md)
  reads a missing coefficient for the tested contrast as a true effect
  of 0, with Type S and Type M undefined, where it stopped with an
  error.

- [`model_data()`](https://pablobernabeu.github.io/pilotr/r/reference/model_data.md)
  and
  [`model_formula()`](https://pablobernabeu.github.io/pilotr/r/reference/model_formula.md)
  read a specification from its path, as documented.

- [`build_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/build_spec.md)
  assembles a valid `exgaussian` response.

- [`build_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/build_spec.md)
  and
  [`spec_from_model()`](https://pablobernabeu.github.io/pilotr/r/reference/spec_from_model.md)
  keep a seed of 2^31 or more, which became `NA`, and the no-code app
  accepts such a seed and prints it in full.

- [`load_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/load_spec.md)
  reads a file that starts with a byte-order mark without a warning.

- [`generate_design_analysis()`](https://pablobernabeu.github.io/pilotr/r/reference/generate_design_analysis.md)
  accepts a focal interaction with its parts in either order and writes
  it under the name brms gives it.

- [`brms_bridge()`](https://pablobernabeu.github.io/pilotr/r/reference/brms_bridge.md)
  and
  [`generate_design_analysis()`](https://pablobernabeu.github.io/pilotr/r/reference/generate_design_analysis.md)
  emit models that brms accepts, with a prior on random-effect standard
  deviations only where the design has random effects and each
  interaction named as brms names it.

- [`run_app()`](https://pablobernabeu.github.io/pilotr/r/reference/run_app.md)
  works without attaching pilotr, keeps a parallel future plan the
  caller has set, restores the caller’s plan when the app stops and
  honours `async = FALSE`.

- The no-code app keeps the session open after an input that its Power
  tab or Verify button cannot use and says what is wrong.

- The no-code app’s summary and plot work after a change of family or
  response name, and the plot draws factors with more than two levels.

- The installed no-code app labels its paste box for screen readers.

### Deprecated

- A within factor whose `vary_within` omits a unit of the design warns,
  since pilotr crosses a within factor with every unit whatever the list
  says. List every unit instead, or make the factor between `item` if
  items carry it. From specification version 0.4, such a list is
  refused.

### Documentation

- The cross-language guarantee is stated as it holds: R and Python
  simulate bit-identical data where they share a maths library, as on
  Linux, except a Beta response whose Gamma shape is exactly 1/3 or 1/4
  under R 4.6.0 or later. It was claimed for the Gaussian family on
  every platform and for any rounded response. The documentation also
  says that one language’s data can differ across operating systems in
  the same way, where it said that a specification and seed always give
  the same data.

- The specification sets out how units are allocated to levels and
  clusters, with worked encodings for a 2 × 2 between design, a two-list
  counterbalanced design and randomisation within clusters, and gives
  every sampling rule in pseudocode, correcting its account of the draws
  the Gamma sampler consumes. Its JSON schema requires each family’s
  parameters and refuses unknown top-level fields and unused response
  parameters. It accepts the single-value forms of `vary_within` and
  `thresholds`, as strict validation does.

- [`?power_curve_mixed`](https://pablobernabeu.github.io/pilotr/r/reference/power_curve_mixed.md),
  [`?solve_curve`](https://pablobernabeu.github.io/pilotr/r/reference/solve_curve.md)
  and the power vignette say that when an effect varies from item to
  item, power at a fixed number of items levels off below 1, however
  many subjects are added, so a curve over subjects can be solved far
  from the true sample size.

- [`?solve_curve`](https://pablobernabeu.github.io/pilotr/r/reference/solve_curve.md)
  reports the coverage of its interval over 1,000 curves a setting,
  where it reported coverage in 35 of 36 checks on curves whose points
  shared their seeds. It says that such a curve, including one built
  from runs that each start from the same seed, gets too narrow an
  interval, and that a finer grid helps only while its points are
  independent.

- The power vignette gains sections on one-sided tests, a preregistered
  simplification rule, the model and test each family gets, a
  sample-size justification with
  [`power_report()`](https://pablobernabeu.github.io/pilotr/r/reference/power_report.md)
  and a comparison of designs with different numbers of subjects and
  items.

- The power and precision vignettes say that their priming example is
  fully crossed and that a counterbalanced study has half the
  observations and is planned with `counterbalance` or the two-list
  encoding. The power vignette and the specification say that simulated
  data are complete, so a study that expects to lose a proportion p of
  its participants recruits N / (1 − p) to analyse N.

- The power and precision vignettes show results computed by this
  version, where their cached results predated 0.3.0, and state that
  each rate is taken over `n_returned`.

- [`?response_variance`](https://pablobernabeu.github.io/pilotr/r/reference/response_variance.md)
  no longer calls the total variance, which includes the fixed effects,
  the right denominator for an effect size, and shows how to turn a
  published Cohen’s d into a coefficient.

- [`?sweep_spec`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md)
  and
  [`?solve_curve`](https://pablobernabeu.github.io/pilotr/r/reference/solve_curve.md)
  say that a sweep over lists,
  [`design_conditions()`](https://pablobernabeu.github.io/pilotr/r/reference/design_conditions.md)
  included, records a grid index, which
  [`solve_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/solve_curve.md)
  solves for unless `x` names another column, where
  [`?sweep_spec`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md)
  said that a sweep could be solved as it stood.

- [`?print.pilotr_power`](https://pablobernabeu.github.io/pilotr/r/reference/print.pilotr_power.md)
  says that `digits` counts decimal places, where it said significant
  digits.

- [`?model_data`](https://pablobernabeu.github.io/pilotr/r/reference/model_data.md)
  and
  [`?precision_design`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md)
  say that
  [`model_data()`](https://pablobernabeu.github.io/pilotr/r/reference/model_data.md)
  logs the response only for the two lognormal families, where
  [`?precision_design`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md)
  said for every family.

- The guides of both builds of the no-code app say that the effect is on
  the scale of the model’s linear predictor, where both said the
  response scale.

- The installed no-code app’s Power tab says that `PILOTR_MAX_SIMS`
  raises its cap on simulations, where it said the number was unlimited
  once the package was installed.

- [`?run_app`](https://pablobernabeu.github.io/pilotr/r/reference/run_app.md)
  and the article on the no-code app say that `async` frees the app’s R
  process for other sessions while the session that started a run waits
  for it, where they said that the app does not block and uses all local
  cores.

## pilotr 0.3.1

CRAN release: 2026-09-20

### Fixes

- [`spec_json()`](https://pablobernabeu.github.io/pilotr/r/reference/spec_json.md)
  no longer corrupts a coefficient on a build of R without long-double
  arithmetic. Numbers were written by asking `jsonlite` for 17
  significant digits and then shortening each one in the text, which
  meant reading the number back with
  [`as.numeric()`](https://rdrr.io/r/base/numeric.html). That read-back
  is not correct everywhere. Without long doubles
  [`as.numeric()`](https://rdrr.io/r/base/numeric.html) accumulates the
  mantissa in a double before applying the decimal exponent, and 17
  significant digits overflow the 53-bit mantissa on the way, so
  `"0.33333333333333331"` came back one unit in the last place below
  `1/3`. The shortened form then recorded that wrong value, and the
  saved specification no longer held the coefficient the user had set.
  Each number is now formatted from the double itself, before `jsonlite`
  sees the document, so nothing is ever read back from the text it was
  just written to. A coefficient of `1/3` is written at 16 significant
  digits rather than 17, which is both the shorter form and the one such
  a build reads back exactly. Found by CRAN’s no-long-double check.

- [`generate_r_script()`](https://pablobernabeu.github.io/pilotr/r/reference/generate_r_script.md)
  could write `Inf` in place of a finite coefficient on the same builds.
  It embeds each number through the helper that looks for the shortest
  decimal that reads back as the same double, and that helper returned
  nothing at all when none of the three widths it tried passed the
  check, a case the caller took to mean the value was infinite. Without
  long doubles the 17-digit width fails that check for the very values
  that need it. The helper now falls back to 17 significant digits,
  which is the widest form there is and so the one to use when nothing
  shorter will do, and returns nothing only for a value that really is
  not finite.

## pilotr 0.3.0

CRAN release: 2026-09-12

### Read this first

Two changes alter numbers that earlier versions produced. Both are
deliberate corrections, and both are unavoidable if the results are to
mean what the documentation says they mean. Install 0.2.1 to reproduce
output from 0.2.1.

- Replicate seeds changed. The power and precision loops previously
  seeded replicate `i` with `seed + (i - 1)`. Consecutive seeds are not
  independent streams in this generator, and the measured consequence
  was severe. The first draw of replicate `i` correlated 0.95 with the
  first draw of replicate `i + 1`. Seeds are now drawn from the shared
  generator, which brings that correlation to -0.02 and moves a
  Ljung-Box test over the replicate means from p below 0.0001 to p of
  0.94. Every number the replicate loops produce is therefore different,
  and the first replicate no longer uses the specification’s own seed.
  [`replicate_seeds()`](https://pablobernabeu.github.io/pilotr/r/reference/replicate_seeds.md)
  is exported so that a hand-written loop or a cluster array task can
  use the same rule.

- An interaction random slope now reaches the data. A slope keyed on an
  interaction, such as `"z_cosine:z_ISI"`, was accepted, sized into the
  covariance, drawn from the random stream and then multiplied by a
  lookup that returned nothing, so it was silently discarded. Meanwhile
  [`model_formula()`](https://pablobernabeu.github.io/pilotr/r/reference/model_formula.md)
  and
  [`brms_bridge()`](https://pablobernabeu.github.io/pilotr/r/reference/brms_bridge.md)
  both emitted that slope, so the model fitted was richer than the
  process that generated the data, which inflates power in the direction
  Barr et al. (2013) describe. Any design that used one was
  anti-conservative and its data have changed.

### Cross-language reproducibility

The R and ‘Python’ implementations agreed on far less than the
documentation claimed. Four separate defects made them diverge, none of
them visible in the shipped examples, because six of the eight set
`response.round`, which quantises exactly these differences away. With
rounding removed, three of the eight examples differed, in up to 9.24%
of cells.

- The accumulation order of the linear predictor is now fixed by the
  specification. R folded the terms one at a time while ‘Python’ summed
  them and added the total. Floating-point addition is not associative,
  so the two disagreed on 63.8% of rows in a matched reproduction over
  200,000 rows.
- Neither language’s built-in summation is used any more. Base R’s
  [`sum()`](https://rdrr.io/r/base/sum.html) accumulates in 80-bit long
  double, and CPython’s [`sum()`](https://rdrr.io/r/base/sum.html) has
  applied Neumaier compensation to floats since 3.12. Both are more
  accurate than a plain double fold, in different ways, so an inner
  product of three terms or more could land on different doubles. Every
  inner product is now an explicit fold.
- Integer powers are written as repeated multiplication. R special-cases
  small integer exponents while ‘Python’ calls the library `pow()`. Over
  200,000 draws in the Gamma sampler’s range the two disagreed on a
  third of inputs by up to 6 ulp. Because that value decides a rejection
  step, it also changed how many draws were consumed. This alone
  accounted for every difference in the `beta_proportion` example.
- The random-effect covariance is bracketed identically. R computed
  `(sd_i · sd_j) · r_ij` and ‘Python’ `(sd_i · r_ij) · sd_j`.
  Multiplication is commutative but not associative, and the difference
  propagated through the Cholesky factor into every random effect drawn.

All eight shipped examples are now bit-identical between the two
languages as shipped, as are five adversarial specifications added to
exercise what the shipped examples do not. With `response.round`
removed, six of the eight remain bit-identical, and the two whose family
applies [`exp()`](https://rdrr.io/r/base/Log.html) to the linear
predictor, `crossed_mixed_rt` and `reading_time_continuous`, differ by
at most one unit in the last place, on 0.09% and 0.05% of cells. That
residue is the libm limit described immediately below, and
`tools/parity/tolerance.json` records which cases carry an allowance for
it, and why.

- The scope of the guarantee is now stated honestly in `spec/SPEC.md`.
  It is exact for `gaussian`, for any design applying no transcendental
  function to the linear predictor, and for any family with
  `response.round` set. The families that apply
  [`exp()`](https://rdrr.io/r/base/Log.html) or
  [`log()`](https://rdrr.io/r/base/Log.html) may differ in the last unit
  in the last place, because IEEE-754 does not require correct rounding
  for those functions and the two builds need not share a maths library.
  Measured rates are given, and the attribution is demonstrated: the
  same design switched to `gaussian`, with an identical seed, structure
  and draw sequence, is bit-identical.

- Three small divergences between the two engines are closed. A
  non-whole seed truncates identically, a whole `spec_version` reads the
  same however it was written, and a non-object unit is now reported as
  such, where R used to crash with a base error. The message text of
  each is character-for-character identical across the twins.

### Validation and versioning

- [`validate_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/validate_spec.md)
  checks a specification against the schema and against the cross-field
  rules the schema cannot express, and
  [`load_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/load_spec.md)
  now calls it by default. A strict draft-07 schema had shipped since
  0.1 with no code path consulting it.

- Validation exists because several ways of getting a specification
  wrong produced plausible data and no error at all. A mistyped
  coefficient key resolved to no column and so silently set that effect
  to zero, which generates exactly the data of a null design and reports
  success. A response parameter left over from another family was
  ignored. Both are now refused.

- `spec_version` negotiates the format. A specification with no such
  field is a 0.2 specification. One using a 0.3 feature must declare
  0.3, because a 0.2 implementation reads it differently and generates
  different data without complaint. One declaring a version newer than
  the implementation understands is refused outright.

- [`simulate_design()`](https://pablobernabeu.github.io/pilotr/r/reference/simulate_design.md)
  gains `validate`, defaulting to `TRUE`. The replicate loops validate
  once and then skip it, so a sweep pays the cost once.

- [`spec_from_model()`](https://pablobernabeu.github.io/pilotr/r/reference/spec_from_model.md)
  gains test coverage: the recovered specification’s units,
  `between`/`vary_within` placement, interaction keys read back off
  product columns, and random-effect estimates are checked against the
  design that generated the pilot data, alongside the refusal paths for
  models the function cannot read.

- The package now declares a minimum R version, `Depends: R (>= 4.0.0)`.
  The cross-language claim assumes the
  [`round()`](https://rdrr.io/r/base/Round.html) that arrived in 4.0.0,
  which is where that floor comes from. A new CI job checks the declared
  Suggests floors, which the matrix (release, devel, oldrel-1) sat well
  above and so never exercised.

- The packaged example specifications are now tested against the
  repository’s own copies. `spec/examples/*.json` is canonical and both
  packages carry a mirror so an installed copy can reach it, but nothing
  enforced the mirror. A load-and-simulate test cannot: a stale packaged
  copy still loads and simulates perfectly well, it simply describes a
  different design from the one the repository documents. The new
  `test-examples.R` compares the bytes, and is twinned with the Python
  suite’s `test_examples.py`. It skips when the package is checked in
  isolation from the repository.

### Fixes

- `n_converged` was not the number of converged fits. The replicate
  loops wrapped each fit in
  [`suppressWarnings()`](https://rdrr.io/r/base/warning.html), so only a
  hard error was visible and every boundary-singular fit counted as
  converged. In the package’s own documented example, 12 subjects by 12
  items with a maximal model, 85 of 100 replicates were
  boundary-singular while `n_converged` reported 100. The loops now
  report `n_attempted`, `n_returned`, `n_converged`, `n_singular` and
  `n_warning` separately. Singular and warning fits are still used,
  since their fixed-effect estimates remain interpretable and dropping
  them would bias the result.

- A non-positive-definite random-effect covariance is now an error,
  naming the grouping factor and the column at which the factorisation
  failed. The failing pivot was previously clamped at zero and the
  result returned in silence, which produced random effects whose
  standard deviations were several times those requested. In one test a
  requested 0.200 came back as 0.805.

- [`spec_json()`](https://pablobernabeu.github.io/pilotr/r/reference/spec_json.md)
  lost precision on nearly every number. It wrote through
  [`as.character()`](https://rdrr.io/r/base/character.html), so a
  coefficient of `1/3` round-tripped to `0.33333333333333298`. Over a
  sample of 214 doubles, 189 failed to round-trip. Since the JSON file
  is the portable artefact, the specification itself was a source of
  divergence. Numbers are now written at the shortest precision that
  round-trips exactly, so `0.3` still reads as `0.3`. The blanket
  `auto_unbox` is gone too, so a one-element `vary_within` or a single
  ordinal threshold stays an array.

- [`generate_r_script()`](https://pablobernabeu.github.io/pilotr/r/reference/generate_r_script.md)
  embedded the specification through
  [`deparse()`](https://rdrr.io/r/base/deparse.html), which prints 15
  significant digits and so does not round-trip. It now emits numbers at
  full precision, which matters because the point of the script is
  bit-for-bit reproduction.

- [`brms_bridge()`](https://pablobernabeu.github.io/pilotr/r/reference/brms_bridge.md)
  wrote its model to the console instead of returning it. The code
  reached standard output on every call, whether or not anything was
  there to read it, and the one internal caller,
  [`generate_design_analysis()`](https://pablobernabeu.github.io/pilotr/r/reference/generate_design_analysis.md),
  had to divert that output around itself through a sink to a temporary
  file. It now returns a `pilotr_bridge` object, visibly, carrying the
  same `formula`, `family`, `priors` and `code` as before, and the new
  [`print()`](https://rdrr.io/r/base/print.html) method writes the code.
  A bare call at the console still shows the model; an assignment is
  silent.

- [`model_formula()`](https://pablobernabeu.github.io/pilotr/r/reference/model_formula.md)
  and
  [`brms_bridge()`](https://pablobernabeu.github.io/pilotr/r/reference/brms_bridge.md)
  emitted correlated random effects unconditionally, along with an LKJ
  prior, while the generative process only correlates them when
  `correlations` is supplied. They now follow the new `correlated` flag
  and emit a double bar otherwise. A group with no slopes keeps a single
  bar, since `lme4` cannot parse `(1 || g)`.

- [`model_data()`](https://pablobernabeu.github.io/pilotr/r/reference/model_data.md)
  did not create the product column an interaction random slope needs,
  so the emitted formula referred to a variable the modelling data
  lacked. It now covers the union of the fixed-coefficient and
  random-slope keys.

- A focal effect that never appears in any fit now warns. It used to
  return decision proportions of zero, which read as ‘this design can
  decide nothing’ when the cause was a name that did not match the
  model.

- A replicate loop in which no fit succeeds now reports why, passing on
  the fitter’s own message. An unidentifiable random-effects structure
  previously produced a silent result of `NA`.

- A correlation naming a random-effect term that does not exist now
  raises a clear error, where it used to fail at a subscript.

- A Poisson mean beyond the sampler’s reach was returned as the
  iteration cap. The inverse-CDF walk starts from `exp(-mean)`, which
  underflows to exactly zero once the mean passes about 746 (a poisson
  intercept of 7 already implies a mean of exp(7), about 1097), after
  which every simulated count came back as the cap of 1000000 while
  reporting success. Both engines now refuse such a mean, naming
  `exp(eta)` and the offending value, with message text
  character-for-character identical across the twins. Feasible means are
  untouched and the parity dumps are unchanged.

- `generate_design_analysis(file = )` wrote the script through a
  text-mode connection. On Windows the file arrived with CRLF line
  endings, which turns the SLURM part’s first line into a
  `#!/bin/bash\r` shebang no cluster can execute, and on every platform
  it carried a doubled trailing newline. The script now reaches the disk
  byte for byte as returned, and the app’s script and specification
  downloads write through the same binary path.

- The emitted SLURM wrapper could only run for its author. It hard-coded
  the author’s cluster account and project paths, so any other user’s
  submission failed at the scheduler while the surrounding instructions
  told them to save the parts under their own names. The wrapper now
  carries two placeholders marked `EDIT`, the `--account` directive and
  a writable `PROJECT_DIR`, and it invokes the `design_analysis.R` that
  sits beside it, wherever the user saved the pair.

- [`precision_design()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md)
  documents all sixteen of its columns. The Monte Carlo standard errors
  and Wilson interval bounds were returned but missing from the
  reference page.

- A true effect of exactly zero returned `Inf`, and the package
  recommends that input.
  [`design_conditions()`](https://pablobernabeu.github.io/pilotr/r/reference/design_conditions.md)
  deliberately produces a null condition so a run can show how often it
  declares something when there is nothing to find. Putting it through
  [`power_design()`](https://pablobernabeu.github.io/pilotr/r/reference/power_design.md)
  returned `type_m = Inf`, straight into the app’s display, while Type S
  silently degenerated to “the estimate is positive” because
  `true_effect > 0` is `FALSE` at zero. The Python twin raised
  `ZeroDivisionError` on the same input. Type S and Type M are now `NA`
  when the true effect is zero or unknown, in both engines, which is the
  guard
  [`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md)
  already carried, applied to the other three sites.

- [`response_variance()`](https://pablobernabeu.github.io/pilotr/r/reference/response_variance.md)
  reported a total that was not the sum of its parts. An undefined fixed
  component was laundered into zero while the result kept calling itself
  the sum. A single-row design now reports a fixed component of 0 and an
  honest total.

- The HPC precision-array job was broken and stale. It hand-listed the
  engine files and omitted `validate.R`, so it died on its first real
  line, and it used the indexed seed rule abandoned at 0.3, which
  correlates consecutive replicates and understates the job’s own Monte
  Carlo error. It now sources the package wholesale, as the parity
  harness does, and uses the package’s own replicate seeds,
  `qnorm(0.975)` and error reporting.

- A test now holds
  [`print.pilotr_power()`](https://pablobernabeu.github.io/pilotr/r/reference/print.pilotr_power.md)
  to keeping its whole output on one stream. A header written with
  [`message()`](https://rdrr.io/r/base/message.html) or through cli
  would land on the message stream while the table beneath it went to
  standard output, which knitr collects separately and renders as two
  boxes for one printed object. The method already used
  [`cat()`](https://rdrr.io/r/base/cat.html) throughout, and
  `test-print.R` keeps it there.

- Every vignette now turns console colour off and fixes the console
  width while it renders. pkgdown passes the calling terminal’s colour
  support into its build subprocess, so a coloured message or error
  would otherwise reach the reader as escape sequences in the middle of
  the text.

### New in the generative core

- `varies_by = "observation"` draws a predictor once per row, for a
  quantity that varies trial by trial. `varies_by` is also validated
  now, because anything other than `"subject"` was previously read as
  item-level, so a predictor declared to vary by `"trial"` was silently
  given one value per item, and anyone who wrote that got wrong results.

- `dist = "uniform"`, with `min` and `max`. A uniform costs the same
  single draw as a normal, so it does not move the stream.

- `reliability` on a predictor simulates imperfect measurement. The
  latent value drives the linear predictor and any slope keyed on it,
  while the contaminated observed value goes into the data. No
  comparable package models unreliable predictors, and cross-level
  interactions are where unreliability bites hardest. Population moments
  are used rather than sample ones, since R’s
  [`mean()`](https://rdrr.io/r/base/mean.html) and
  [`sd()`](https://rdrr.io/r/stats/sd.html) accumulate in long double
  and ’Python’s do not.

- The `exgaussian` family, the registered model family for reaction-time
  work, in brms’s parameterisation so that a specification and the model
  fitted to it agree on what the intercept means. A shifted lognormal is
  not a substitute, because
  [`model_data()`](https://pablobernabeu.github.io/pilotr/r/reference/model_data.md)
  logs the response back and leaves a symmetric residual on the analysis
  scale.

- [`response_variance()`](https://pablobernabeu.github.io/pilotr/r/reference/response_variance.md)
  decomposes the linear predictor’s variance into the fixed part, each
  grouping factor’s part, and the residual. Each grouping factor’s
  component is exact for the realised design, averaging over the
  random-effect distribution analytically, because estimating it from
  the drawn effects of 30 subjects carries a sampling error of around a
  quarter of the component.

  A residual is reported for all eight families, not only the four
  carrying an explicit `sigma`, so the components are a complete
  decomposition everywhere and their ratios read as the design’s
  intraclass correlations. For the link families it is the latent-scale
  distribution-specific variance (Nakagawa, Johnson and Schielzeth,
  2017). Three of those four are derived from the process pilotr
  simulates and are exact for it: a `bernoulli` row is drawn as
  `1[u < invlogit(eta)]`, so its latent error is a standard logistic
  variate of variance `pi^2 / 3`, `ordinal` compares the same uniform
  against cumulative thresholds and inherits it, and for `beta` the
  identity `Var(logit(Y)) = trigamma(a) + trigamma(b)` is exact.
  `poisson` is the one approximation, since a count of zero has no
  logarithm. The trigamma form is used, and it is worth reading as an
  order of magnitude when counts are rare, where the published
  alternatives diverge sharply from it.

- [`calibrate_response()`](https://pablobernabeu.github.io/pilotr/r/reference/calibrate_response.md)
  rescales a design to a target total variance, which is what lets a
  region of practical equivalence be stated in standard-deviation units
  and read the same way across designs. It now accounts for a residual
  it cannot move. A `bernoulli` or `ordinal` design carries a latent
  residual of about 3.29, so calibrating one to a total variance of 1 is
  impossible, and is refused with that number in the message. For
  `poisson` and `beta` the residual moves with the linear predictor, so
  the factor is solved numerically, and the whole search costs one extra
  simulation.

### New in the design-analysis layer

- [`generate_design_analysis()`](https://pablobernabeu.github.io/pilotr/r/reference/generate_design_analysis.md)
  emits a Bayesian design analysis as a runnable script: a
  [`brm()`](https://paulbuerkner.com/brms/reference/brm.html) call with
  `sample_prior = "yes"`, a Savage-Dickey Bayes factor, a
  highest-density interval against a region of practical equivalence, a
  three-way supported/null/inconclusive verdict, and a convergence gate
  that withholds every verdict when R-hat or the divergence rate fails
  it. Optionally it also emits a SLURM array wrapper and an aggregator.
  Emitting a script is what keeps this reachable from the browser build,
  where Stan cannot run at all.

- [`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md)
  is no longer restricted to one within-unit factor and a crossed
  design, and no longer fits a formula written into the source. It takes
  `focal` and `formula`, derives the model from the specification
  through
  [`model_formula()`](https://pablobernabeu.github.io/pilotr/r/reference/model_formula.md)
  and
  [`model_data()`](https://pablobernabeu.github.io/pilotr/r/reference/model_data.md),
  and runs the same replicate loop as
  [`precision_design()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md).
  A design with two factors, with continuous predictors, or with a
  smaller random-effects structure than the maximal one was previously
  refused outright or analysed under a model it had not described.

- [`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md)
  runs an analysis over any addressed field of a specification.
  [`power_curve_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_curve_mixed.md)
  and
  [`precision_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_curve.md)
  are now thin wrappers over it for sample size, which was previously
  the only axis reachable without writing the loop by hand.
  [`design_conditions()`](https://pablobernabeu.github.io/pilotr/r/reference/design_conditions.md)
  builds the coefficient sets for an effect-size sweep, including a
  shared all-zero condition.

- Every reported rate now carries its Monte Carlo standard error and a
  Wilson interval. The curve functions default to 60 replicates, at
  which a rate near 0.5 has a standard error of 0.065, which cannot
  support a claim anyone would want to make. The Wilson interval is used
  because the plain standard error collapses to zero at rates of 0 and
  1, reading as certainty precisely where there is least of it.

- [`precision_design()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md)
  and
  [`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md)
  default `focal` to every coefficient in the specification.

- [`print()`](https://rdrr.io/r/base/print.html) for a power result
  shows the fit accounting and each estimate beside its uncertainty.

- [`solve_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/solve_curve.md)
  and
  [`target_n()`](https://pablobernabeu.github.io/pilotr/r/reference/target_n.md)
  solve a simulated design curve for the value that meets a target. The
  package computed the whole curve and then handed back the last and
  most consequential step, the number that goes into a preregistration,
  at replicate counts where neighbouring points are not significantly
  different. The plots drew a dashed line at 0.80 and left the reader to
  judge the crossing by eye.
  [`target_n()`](https://pablobernabeu.github.io/pilotr/r/reference/target_n.md)
  now takes the data frame
  [`power_curve_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_curve_mixed.md)
  already returns and reports the sample size at which power reaches the
  target, with a confidence interval on it, rounded up to whole
  subjects.
  [`solve_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/solve_curve.md)
  is the general form and reads
  [`precision_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_curve.md)
  and
  [`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md)
  output through the same column names, on any axis the sweep varied.

  The fit is a binomial probit regression of the decision rate against
  the swept value, weighted by the replicates behind each point and
  inverted by the delta method that
  [`MASS::dose.p()`](https://rdrr.io/pkg/MASS/man/dose.p.html) applies
  to a fitted `glm`. The probit follows from the design itself: under
  the normal approximation to a two-group comparison the probit of power
  is linear in the square root of the sample size. Checked against
  [`stats::power.t.test()`](https://rdrr.io/r/stats/power.t.test.html)
  over 36 combinations of effect size, target power and grid shape, at
  400 replicates a point, the solved sample size sat within 2.9% of the
  analytic answer on average and 6.9% at worst, against 3.4% and 8.9%
  for a logit fitted the same way.
  `tools/calibration/solve_curve_calibration.R` runs that comparison and
  writes every figure quoted for it to a file beside itself. Where the
  two-parameter model does not describe the curve, the interval is
  widened by the heterogeneity factor of probit analysis, reported as
  `dispersion`, since a curve the model fits badly should not report the
  narrow interval its replicate counts alone imply.

  Nothing extrapolates. A curve that does not reach the target within
  the sizes it swept is refused, with the range it did cover, and so is
  a fit that solves past the end of the sweep, a curve with no trend to
  invert, and a slope that cannot be told from zero. The app reports the
  solved size beside the curve it draws, and prints the refusal when the
  sweep does not settle the question.

- The Bayesian design-analysis record carries the rule that produced it,
  namely the decision thresholds, the convergence gate, a specification
  fingerprint and the package version, and the aggregator refuses to
  pool replicates that disagree on any of them.

### References

Barr, D. J., Levy, R., Scheepers, C., & Tily, H. J. (2013). Random
effects structure for confirmatory hypothesis testing: Keep it maximal.
*Journal of Memory and Language, 68*(3), 255-278.
<doi:10.1016/j.jml.2012.11.001>

Nakagawa, S., Johnson, P. C. D. and Schielzeth, H. (2017). The
coefficient of determination R2 and intra-class correlation coefficient
from generalized linear mixed-effects models revisited and expanded.
*Journal of the Royal Society Interface, 14*(134), 20170213.
<doi:10.1098/rsif.2017.0213>

## pilotr 0.2.1

### Documentation

- A package-level help topic,
  [`?pilotr`](https://pablobernabeu.github.io/pilotr/r/reference/pilotr-package.md),
  introduces the workflow and groups every export by the stage it
  belongs to.
- The reference examples that cut the replicate count so they run
  quickly now say so, and each `n_sims` argument documents the Monte
  Carlo error and a usable recommendation for real planning.
- The
  [`precision_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_curve.md)
  example passes a region of practical equivalence clearly inside the
  effect it declares. It previously fell back to a default equal to that
  effect, the configuration the vignette warns against.

### Fixes

- [`model_formula()`](https://pablobernabeu.github.io/pilotr/r/reference/model_formula.md)
  binds the global environment, so a formula prints without the
  environment tag that varied between builds and leaked into the
  rendered reference pages.

## pilotr 0.2.0

### Design specifications

- [`pilotr_example()`](https://pablobernabeu.github.io/pilotr/r/reference/pilotr_example.md)
  lists the design specifications shipped with the package, one per
  design family, and returns the path to each for
  [`load_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/load_spec.md).
  They are the same JSON files that drive the Python package and the
  no-code app, and the Python package gains the same function.
- The eight per-family specifications now travel inside the package
  itself, so an installed copy can load them without a checkout of the
  repository.

### Documentation

- A new *Worked examples* article simulates every shipped specification.

## pilotr 0.1.0

Initial release.

### Design specifications and simulation

- [`simulate_design()`](https://pablobernabeu.github.io/pilotr/r/reference/simulate_design.md)
  generates a data set from a portable JSON design specification shared
  with the Python package of the same name, with
  [`build_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/build_spec.md)
  composing a specification from a flat list of design inputs,
  [`load_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/load_spec.md)
  reading one authored elsewhere (such as one downloaded from the
  no-code app) and
  [`spec_json()`](https://pablobernabeu.github.io/pilotr/r/reference/spec_json.md)
  serialising it back to JSON.
- The response families are Gaussian, lognormal, shifted lognormal,
  Bernoulli, Poisson, ordinal and Beta.
- A shared cross-language random-number generator
  ([`make_rng()`](https://pablobernabeu.github.io/pilotr/r/reference/make_rng.md),
  with the AS 241 inverse normal in
  [`as241()`](https://pablobernabeu.github.io/pilotr/r/reference/as241.md))
  makes the simulated data bit-identical to the Python package’s given
  the same specification and seed.
- A `per_subject` value below 1 or above the number of items is rejected
  with a clear error. The Python package raises on the same inputs.

### Power and precision analysis

- Simulation-based power and design analysis with
  [`power_design()`](https://pablobernabeu.github.io/pilotr/r/reference/power_design.md),
  [`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md)
  and
  [`power_curve_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_curve_mixed.md),
  reporting power alongside the Type S and Type M errors of Gelman and
  Carlin (2014).
  [`power_curve_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_curve_mixed.md)
  returns an `n_converged` column alongside `power` and `type_m`, so the
  Monte Carlo standard error can be computed over the converged
  replicates.
- Precision and region-of-practical-equivalence (ROPE) design analysis
  with
  [`precision_design()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md)
  and
  [`precision_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_curve.md).
- Every analysis seeds the shared RNG as `seed + (replicate - 1)`, the
  same indexed-seed rule across all power and precision functions, the
  HPC array script and both languages.
- [`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md)
  and
  [`power_curve_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_curve_mixed.md)
  reject a specification without an item unit with a clear error instead
  of returning `NaN` power.
- All simulation-based power and precision analyses take a `workers`
  argument that spreads the Monte Carlo replicates across local cores
  with base R’s `parallel` package. Because every replicate seeds the
  shared RNG from its own index, any worker count returns results
  identical to a serial run, and the sweep functions create their worker
  pool once and reuse it across all sample sizes.

### Modelling bridges and the no-code app

- [`model_formula()`](https://pablobernabeu.github.io/pilotr/r/reference/model_formula.md)
  and
  [`model_data()`](https://pablobernabeu.github.io/pilotr/r/reference/model_data.md)
  derive the `lmer` formula and modelling data frame implied by a
  specification, and
  [`brms_bridge()`](https://pablobernabeu.github.io/pilotr/r/reference/brms_bridge.md)
  derives a `brms` formula, family and priors, mapping every response
  family to its `brms` counterpart (including
  [`Beta()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)).
- [`generate_r_script()`](https://pablobernabeu.github.io/pilotr/r/reference/generate_r_script.md)
  writes a self-contained, reproducible R script from a specification.
- A no-code Shiny application over the same specification, launched with
  [`run_app()`](https://pablobernabeu.github.io/pilotr/r/reference/run_app.md).

### Metadata

- `citation("pilotr")` builds its version note from the package
  metadata, so it cannot drift from `DESCRIPTION`.
