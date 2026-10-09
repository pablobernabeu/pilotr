# Changelog

## pilotr (development version)

### Read this first

- [`power_design()`](https://pablobernabeu.github.io/pilotr/r/reference/power_design.md)
  refuses designs whose rows are correlated. It t-tests every row as an
  independent observation, yet it accepted any Gaussian design with a
  two-level between factor. That included designs crossed with items,
  designs with a within factor and subjects nested in sites. Take 30
  subjects crossed with 20 items, with by-subject and residual standard
  deviations of 1 and a by-item one of 0.3. With no true effect, the
  row-wise t-test was significant in 105 of 200 replicates, against 16
  for a t-test of the subject means on the same data. The power it
  reported was inflated and its Type M understated. It now takes one row
  per subject and no clustering, and otherwise names the item unit,
  within factor or grouping factor it found and points to
  [`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md).
  The no-code app, installed or in the browser, applies the same rule.

- [`brms_bridge()`](https://pablobernabeu.github.io/pilotr/r/reference/brms_bridge.md)
  and
  [`generate_design_analysis()`](https://pablobernabeu.github.io/pilotr/r/reference/generate_design_analysis.md)
  emit models that brms accepts, with priors on the data’s scale. A
  prior on random-effect standard deviations was emitted for designs
  without random effects, which brms rejects, so four of the eight
  shipped examples failed before sampling. An interaction key was used
  as written, `cond:age`, where brms names the coefficient `age:cond`,
  so brms refused a fifth example’s prior on it. The priors were also
  fixed at a unit scale described as standardised, although nothing
  standardised the data. On a Gaussian design with a residual SD of 10,
  the effect’s prior was N(0, 0.5) and the Savage-Dickey Bayes factor
  was about 1 whatever the data. Priors for the continuous families now
  scale with the response and with each column’s realised standard
  deviation, the intercept takes brms’s default, and the emitted call
  says that its data come from
  [`model_data()`](https://pablobernabeu.github.io/pilotr/r/reference/model_data.md).
  Emitted scripts, priors and Bayes factors change.

- [`generate_design_analysis()`](https://pablobernabeu.github.io/pilotr/r/reference/generate_design_analysis.md)
  accepts a focal interaction with its parts in either order and writes
  it into the script under the name brms gives it.
  [`brms_bridge()`](https://pablobernabeu.github.io/pilotr/r/reference/brms_bridge.md)
  validates a specification given as a list, as the other entry points
  do, since its priors are now measured on the data that specification
  simulates. It refuses a `prior_scale` or `interaction_scale` that is
  not a single positive number, and a specification that gives one
  interaction two keys, such as `cond:age` and `age:cond`, since brms
  estimates one coefficient for both. A coefficient whose column does
  not vary is reported in a warning.

- A specification with two or more factors between the same unit is
  refused. pilotr assigns the levels of each between factor to blocks of
  units on its own, so two such factors fell into the same or
  overlapping blocks. A 2 × 2 between-subjects design over 40 subjects
  produced cells of 20, 0, 0 and 20, and the second effect and the
  interaction could not be estimated. The refusal names the encoding
  that works in every version, one between factor whose levels are the
  cells.
  [`spec_from_model()`](https://pablobernabeu.github.io/pilotr/r/reference/spec_from_model.md)
  now builds that factor itself when it reads two factors between one
  unit off a fitted pilot, and says so in a message. It used to write a
  specification that left half the cells of a 2 × 2 pilot empty.

- [`spec_from_model()`](https://pablobernabeu.github.io/pilotr/r/reference/spec_from_model.md)
  warns when the pilot was counterbalanced, each subject seeing each
  item under one condition, because the specification it returns crosses
  every subject with every item under every condition. Planning a
  counterbalanced study that way doubles the observations per subject
  and overstates power. By the formulas of Judd, Westfall and Kenny
  (2017), the `crossed_mixed_rt` example, with 30 subjects, 24 items and
  an effect of 0.05 on the log scale, has a power of about 0.80 crossed
  and 0.55 counterbalanced. The warning names the encoding that plans a
  counterbalanced study, a list factor between subjects and an item-set
  factor between items. The specification returned is unchanged.

- [`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md)
  and
  [`precision_design()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md)
  fit the model each family needs. They fitted a linear mixed model to
  every response, so for `bernoulli` and `poisson` outcomes the estimate
  was a difference on the response scale, compared with a true value on
  the link scale. Take 20 subjects crossed with 16 items and an effect
  of 0.5 on the logit scale. Over 30 replicates, the mean estimate was
  0.090 and Type M 0.21, where a logistic mixed model of the same
  replicates gives 0.47 and 1.12, so real exaggeration read as
  underestimation. Against a region of practical equivalence of 0.1, the
  interval was 0.14 wide and fell outside the region in none of the
  replicates, against 0.71 wide and 17 of 30 on the logit scale. Counts
  on the same crossing, with an effect of 0.3 on the log scale, gave a
  mean estimate of 0.90 and Type M 3.01 over 20 replicates, against 0.31
  and 1.02. Bernoulli and Poisson designs are now fitted with
  [`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) and tested
  with a Wald z, so Type M, mean estimates, interval widths and ROPE
  decisions change materially. Power barely moves: 22 of the 30 accuracy
  replicates are significant under either model. A Wald z is the test
  Bolker et al. (2009) recommend for such a model. With few subjects or
  items, it rejects more often than its nominal level (Li and Redden,
  2015).

- A design without random effects, which failed every replicate with “No
  random effects terms specified in formula”, is fitted with
  [`lm()`](https://rdrr.io/r/stats/lm.html) or
  [`glm()`](https://rdrr.io/r/stats/glm.html), so every between-subjects
  example has a power path. That includes the shipped
  `poisson_counts_between`, `ordinal_likert_between` and
  `beta_proportion`, which had none in either language.

- For `ordinal` and `beta` outcomes, which have no frequentist backend
  here yet, both functions warn, withhold the mean estimate, Type M and
  the ROPE decisions, and point to
  [`generate_design_analysis()`](https://pablobernabeu.github.io/pilotr/r/reference/generate_design_analysis.md).
  A later release will refuse them or fit a proportional-odds and a Beta
  model. Type S, which compares signs, is kept for main effects and
  withheld for interactions. A bounded response compresses differences
  near the ends of its range, so an interaction that is positive on the
  logit scale can be negative on the response scale.

- Results name the fitter and the test.
  [`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md)
  returns them as `fitter` and prints them, and
  [`precision_design()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md)
  and the sweeps carry them as a column. A formula given by the user is
  fitted as written, by `lmer()`, or by
  [`lm()`](https://rdrr.io/r/stats/lm.html) when it has no random terms,
  whatever the family.

### Fixes

- The no-code app, installed or in the browser, sends a design its power
  tab cannot analyse to the analysis that fits it. It sent every such
  design to the
  [`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md)
  of both packages. The Python one is now offered only for one within
  factor crossed with items, in a family whose coefficients are on the
  scale it analyses, and ordinal and Beta designs are sent to
  [`generate_design_analysis()`](https://pablobernabeu.github.io/pilotr/r/reference/generate_design_analysis.md)
  and
  [`brms_bridge()`](https://pablobernabeu.github.io/pilotr/r/reference/brms_bridge.md).

- `n_warning` no longer counts boundary-singular fits, which
  `n_singular` already counts. lme4 records its singular-fit notice
  among the optimiser’s messages, so every singular fit was also counted
  as a fit with warnings, and a run could not show whether any fit had
  failed to converge. `n_converged` and `n_singular` are unchanged.

- Every specification field is read by its exact name. R looked fields
  up with `$`, which matches a name partially, so a design whose random
  entries were `subject_site` and no `subject` acquired by-subject
  random effects it never declared, and its data differed from the
  Python twin’s in every row. With `strict = FALSE`, a misspelt field
  such as `rounding` acted as the field it resembled, although the
  warning said it would be ignored.

- [`load_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/load_spec.md)
  and
  [`validate_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/validate_spec.md)
  refuse a specification that repeats a key within one object. JSON
  leaves a repeated key undefined and the two readers disagree:
  `{"grp": 0.3, "grp": 0.5}` applied 0.3 twice in R and 0.5 once in
  Python. The refusal does not depend on `validate`.

- A seed must lie within ±(2^53 − 1), the range in which both twins read
  a JSON integer exactly. Larger seeds were read as different numbers by
  the two readers, which then drew different data.

- [`build_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/build_spec.md)
  and
  [`spec_from_model()`](https://pablobernabeu.github.io/pilotr/r/reference/spec_from_model.md)
  keep a seed of 2^31 or more, which they used to turn into `NA`, and
  leave a non-whole seed for
  [`validate_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/validate_spec.md)
  to refuse, where they used to truncate it. Seeds in the integer range
  are stored as before. The app accepts a larger seed for the same
  reason and prints it in full.

- A factor level, a `vary_within` entry or a `correlated` flag written
  as `null` is refused, as the Python twin refuses it. R labelled those
  rows `NA` and read the flag as `false`.

- A malformed `response`, `predictors`, `random` or `spec_version` value
  is reported in the validator’s own words. Some shapes stopped with a
  base R error, and a `spec_version` given as an array of two versions
  was read as the first.

- [`load_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/load_spec.md)
  reads a file that starts with a byte-order mark without warning about
  it. The specification format allows one, since some Windows editors
  add it.

- [`design_conditions()`](https://pablobernabeu.github.io/pilotr/r/reference/design_conditions.md)
  changes only the effects it names. Through
  [`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md),
  its conditions replaced the whole coefficient list, so every
  coefficient they did not name was dropped from data generation, from
  [`model_formula()`](https://pablobernabeu.github.io/pilotr/r/reference/model_formula.md)
  and from the default focal set, contrary to the documentation. On
  `reading_time_continuous`, a sweep of `cond` simulated and analysed a
  one-coefficient model. Conditions are now merged into the
  specification’s own coefficients, in its order, and a condition naming
  a coefficient the specification lacks is refused. They carry the class
  `pilotr_conditions`, which subsetting them with `[` and joining them
  with [`c()`](https://rdrr.io/r/base/c.html) keep. Any other list still
  replaces the addressed field wholesale. With 0.3.x, pass
  `.base = spec$fixed$coefficients` to
  [`design_conditions()`](https://pablobernabeu.github.io/pilotr/r/reference/design_conditions.md)
  for the same result.

- The
  [`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md)
  example no longer runs the null condition twice.

- [`pilotr::run_app()`](https://pablobernabeu.github.io/pilotr/r/reference/run_app.md)
  works in a session that has not attached pilotr. The app looked for
  the package’s functions among attached packages and then for source
  files, and in an installed package found neither, so the documented
  launcher stopped with “no pilotr engine sources”. The app now takes
  the functions from the package’s namespace, without attaching the
  package.

- [`validate_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/validate_spec.md)
  refuses specifications whose names collide. Each of these validated
  and then moved an effect, rescaled a variance or overwrote a column
  without a word: two columns with one name (a factor called `subject`,
  a response named like a factor), a contrast column defined by two
  factors or named like a predictor, a unit, a grouping factor, another
  factor or the response, an interaction whose analysis column `a_b`
  already exists, a level listed twice, a correlation pairing a term
  with itself or giving one pair twice, a factor both between and within
  a unit, and a `random.item` entry in a design without items, which was
  dropped. Each refusal reads the same in the Python twin. The rules are
  set out under Names in the specification.

- [`validate_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/validate_spec.md)
  refuses a factor, predictor or extra grouping-factor name that is
  empty or made only of spaces, and a response name made only of spaces
  (an empty one was already refused). Such a name validated and then
  became a column with no visible name. An empty factor name stopped
  [`simulate_design()`](https://pablobernabeu.github.io/pilotr/r/reference/simulate_design.md)
  with base R’s “replacement has length zero”, which does not say which
  field is at fault, and a blank `random` key was reported as
  `random. must be an object`.

- Two refusals read as the Python twin’s do. A `correlated` flag that is
  not a boolean must be `true or false`, as the JSON file spells it, and
  `response.thresholds` must be “a number or a non-empty numeric array”.

- [`model_data()`](https://pablobernabeu.github.io/pilotr/r/reference/model_data.md)
  no longer fills a factor’s later contrast columns with `NA` when the
  first contrast column carries the factor’s own name. With three
  levels, every
  [`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md)
  replicate then failed to fit.

- [`power_design()`](https://pablobernabeu.github.io/pilotr/r/reference/power_design.md)
  reports a true effect of 0, with Type S and Type M undefined (`NaN`),
  when the specification has no coefficient for the factor’s contrast,
  as with `"coefficients": {}`. It stopped with “missing value where
  TRUE/FALSE needed”.

### Deprecated

- A within factor whose `vary_within` omits a unit of the design now
  warns, in strict and lenient validation alike. pilotr has always
  crossed a within factor with every unit, so the list changed nothing.
  From spec version 0.4 such a list will be refused.

### Documentation

- `citation("pilotr")` names the R package version alone. It read “R and
  Python package version 0.3.1” while the Python package stood at 0.3.0,
  and the two packages are released separately. It also gives the
  author’s ORCID iD now.

- The 0.3.0 and 0.3.1 entries below no longer give release dates. Both
  gave dates on which neither version reached CRAN, and the site already
  prints the CRAN date beside each entry.

- The tree CRAN published as pilotr 0.3.0 is tagged `cran-0.3.0` on
  GitHub, and 0.3.1 is tagged `v0.3.1` and archived on Zenodo. The
  `v0.3.0` tag and its Zenodo record predate two changes that CRAN’s
  0.3.0 carried, namely that
  [`brms_bridge()`](https://pablobernabeu.github.io/pilotr/r/reference/brms_bridge.md)
  returns a visible `pilotr_bridge` object and that this object has a
  print method.

- [`sweep_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/sweep_spec.md)
  and
  [`solve_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/solve_curve.md)
  say that a sweep over lists,
  [`design_conditions()`](https://pablobernabeu.github.io/pilotr/r/reference/design_conditions.md)
  included, records a grid index in its leading column, and that
  [`solve_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/solve_curve.md)
  solves for that index unless `x` names another column. Both pages show
  the sweep to use for a single effect,
  `sweep_spec(spec, "fixed$coefficients$cond", values, power_mixed)`. It
  keeps the other coefficients, and
  [`solve_curve()`](https://pablobernabeu.github.io/pilotr/r/reference/solve_curve.md)
  returns a minimum detectable effect from its curve.

- The specification states how units are allocated: the block rule for
  between factors and extra grouping factors, the nesting of a between
  factor within clusters that follows from it, and the full crossing of
  within factors. It gives worked encodings for a 2 × 2 between design,
  a two-list counterbalanced design and randomisation within clusters,
  each held bit-identical across the twins by the parity harness. The
  power vignette says that its priming example is fully crossed, and
  that simulated data are complete, so expected attrition is allowed for
  by recruiting N / (1 − p).

- The precision vignette says that its priming example is fully crossed,
  and that a counterbalanced study has half the observations and is
  planned with the two-list encoding.
  [`?spec_from_model`](https://pablobernabeu.github.io/pilotr/r/reference/spec_from_model.md)
  explains the counterbalancing warning and that encoding. It also shows
  a `prep` function for
  [`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md)
  and
  [`precision_design()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md)
  that drops the rows a counterbalanced study would not observe and runs
  in the worker processes of `workers > 1`.

### References

Bolker, B. M., Brooks, M. E., Clark, C. J., Geange, S. W., Poulsen, J.
R., Stevens, M. H. H. and White, J.-S. S. (2009). Generalized linear
mixed models: A practical guide for ecology and evolution. *Trends in
Ecology & Evolution, 24*(3), 127-135. <doi:10.1016/j.tree.2008.10.008>

Judd, C. M., Westfall, J. and Kenny, D. A. (2017). Experiments with more
than one random factor: Designs, analytic models, and statistical power.
*Annual Review of Psychology, 68*, 601-625.
<doi:10.1146/annurev-psych-122414-033702>

Li, P. and Redden, D. T. (2015). Comparing denominator degrees of
freedom approximations for the generalized linear mixed model in
analyzing binary outcome in small sample cluster-randomized trials. *BMC
Medical Research Methodology, 15*, 38. <doi:10.1186/s12874-015-0026-x>

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
  emits a Bayesian design analysis as a runnable script: a `brm()` call
  with `sample_prior = "yes"`, a Savage-Dickey Bayes factor, a
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

Released 2026-07-23.

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

Released 2026-07-15.

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
  family to its `brms` counterpart (including `Beta()`).
- [`generate_r_script()`](https://pablobernabeu.github.io/pilotr/r/reference/generate_r_script.md)
  writes a self-contained, reproducible R script from a specification.
- A no-code Shiny application over the same specification, launched with
  [`run_app()`](https://pablobernabeu.github.io/pilotr/r/reference/run_app.md).

### Metadata

- `citation("pilotr")` builds its version note from the package
  metadata, so it cannot drift from `DESCRIPTION`.
