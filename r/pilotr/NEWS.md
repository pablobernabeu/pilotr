# pilotr (development version)

## Read this first

* `power_design()` refuses designs whose rows are correlated. It t-tests every row as an
  independent observation, yet it accepted any Gaussian design with a two-level between factor.
  That included designs crossed with items, designs with a within factor and subjects nested in
  sites. Take 30 subjects crossed with 20 items, with by-subject and residual standard
  deviations of 1 and a by-item one of 0.3. With no true effect, the row-wise t-test was
  significant in 105 of 200 replicates, against 16 for a t-test of the subject means on the same
  data. The power it reported was inflated and its Type M understated. It now takes one row per
  subject and no clustering, and otherwise names the item unit, within factor or grouping factor
  it found and points to `power_mixed()`. The no-code app, installed or in the browser, applies
  the same rule.

* `brms_bridge()` and `generate_design_analysis()` emit models that brms accepts, with priors on
  the data's scale. A prior on random-effect standard deviations was emitted for designs without
  random effects, which brms rejects, so four of the eight shipped examples failed before
  sampling. An interaction key was used as written, `cond:age`, where brms names the coefficient
  `age:cond`, so brms refused a fifth example's prior on it. The priors were also fixed at a unit
  scale described as standardised, although nothing standardised the data. On a Gaussian design
  with a residual SD of 10, the effect's prior was N(0, 0.5) and the Savage-Dickey Bayes factor
  was about 1 whatever the data. Priors for the continuous families now scale with the response
  and with each column's realised standard deviation, the intercept takes brms's default, and
  the emitted call says that its data come from `model_data()`. Emitted scripts, priors and
  Bayes factors change.

* `generate_design_analysis()` accepts a focal interaction with its parts in either order and
  writes it into the script under the name brms gives it. `brms_bridge()` validates a
  specification given as a list, as the other entry points do, since its priors are now
  measured on the data that specification simulates. It refuses a `prior_scale` or
  `interaction_scale` that is not a single positive number, and a specification that gives one
  interaction two keys, such as `cond:age` and `age:cond`, since brms estimates one coefficient
  for both. A coefficient whose column does not vary is reported in a warning.

* A specification with two or more factors between the same unit is refused unless it declares
  specification version 0.4, which allocates such factors jointly. pilotr assigned the levels of
  each between factor to blocks of units on its own, so two such factors fell into the same or
  overlapping blocks. A 2 × 2 between-subjects design over 40 subjects produced cells of 20, 0, 0
  and 20, and the second effect and the interaction could not be estimated. The refusal names the
  encoding that works in every version, one between factor whose levels are the cells.
  `spec_from_model()` now builds that factor itself when it reads two factors between one unit
  off a fitted pilot, and says so in a message. It used to write a specification that left half
  the cells of a 2 × 2 pilot empty.

* Specification version 0.4 allocates several between factors on one unit jointly, so a 2 × 2
  between-subjects design over 40 subjects has four cells of 10. The units fall into one block
  per combination of levels, with the first-listed factor varying slowest, and nothing is drawn.
  A specification using this must declare `"spec_version": "0.4"`, a version pilotr 0.3.x
  refuses to read. One that does not still gets the refusal introduced above, which names both
  remedies. `build_spec()` and `spec_from_model()` declare the lowest version a specification
  needs, so their output stays readable by earlier releases unless it uses a 0.4 feature. A 0.4
  specification whose `vary_within` omits a unit of the design is refused.

* `spec_from_model()` turned a counterbalanced pilot, each subject seeing each item under one
  condition, into a specification that crosses every subject with every item under every
  condition, without saying so. Planning a counterbalanced study that way doubles the observations
  per subject and overstates power. By the formulas of Judd, Westfall and Kenny (2017), the
  `crossed_mixed_rt` example, with 30 subjects, 24 items and an effect of 0.05 on the log scale,
  has a power of about 0.80 crossed and 0.55 counterbalanced. The specification it returns for
  such a pilot now sets `counterbalance` on the factor, a field of specification version 0.4 (see
  New features). `counterbalance = FALSE` returns the crossed specification, which earlier
  releases read, with a warning. The warning names the field, and for a two-level factor also
  the encoding that earlier releases read, a list factor between subjects and an item-set factor
  between items.

* `power_mixed()` and `precision_design()` fit the model each family needs. They fitted a linear
  mixed model to every response, so for `bernoulli` and `poisson` outcomes the estimate was a
  difference on the response scale, compared with a true value on the link scale. Take 20
  subjects crossed with 16 items and an effect of 0.5 on the logit scale. Over 30 replicates, the
  mean estimate was 0.090 and Type M 0.21, where a logistic mixed model of the same replicates
  gives 0.47 and 1.12, so real exaggeration read as underestimation. Against a region of
  practical equivalence of 0.1, the interval was 0.14 wide and fell outside the region in none of
  the replicates, against 0.71 wide and 17 of 30 on the logit scale. Counts on the same crossing,
  with an effect of 0.3 on the log scale, gave a mean estimate of 0.90 and Type M 3.01 over 20
  replicates, against 0.31 and 1.02. Bernoulli and Poisson designs are now fitted with
  `lme4::glmer()` and tested with a Wald z, so Type M, mean estimates, interval widths and ROPE
  decisions change materially. Power barely moves: 22 of the 30 accuracy replicates are
  significant under either model. A Wald z is the test Bolker et al. (2009) recommend for such a
  model. With few subjects or items, it rejects more often than its nominal level (Li and Redden,
  2015).

* A design without random effects, which failed every replicate with "No random effects terms
  specified in formula", is fitted with `lm()` or `glm()`, so every between-subjects example has
  a power path. That includes the shipped `poisson_counts_between`, `ordinal_likert_between` and
  `beta_proportion`, which had none in either language.

* For `ordinal` and `beta` outcomes, which have no frequentist backend here yet, both functions
  warn, withhold the mean estimate, Type M and the ROPE decisions, and point to
  `generate_design_analysis()`. A later release will refuse them or fit a proportional-odds and a
  Beta model. Type S, which compares signs, is kept for main effects and withheld for
  interactions. A bounded response compresses differences near the ends of its range, so an
  interaction that is positive on the logit scale can be negative on the response scale.

* Results name the fitter and the test. `power_mixed()` returns them as `fitter` and prints them,
  and `precision_design()` and the sweeps carry them as a column. A formula given by the user is
  fitted as written, by `lmer()`, or by `lm()` when it has no random terms, whatever the family.

* `calibrate_response(tune = "all")` holds the intercept and rescales only what carries
  variance, since the intercept is a location. Multiplying the intercept too raised every
  response of the lognormal families, less any shift, to the power of the scaling factor. The
  documented example turned a median reaction time of 600 ms into 16.6 hours. For `poisson` and
  `beta`, the search for the factor began where every term, the intercept included, is zero.
  `poisson_counts_between`, whose own total variance is 0.27, was refused a target of 1 by a
  message that put the family's residual alone at 1.64. A target of 2 returned an intercept of
  11.4, too large to simulate. The example now calibrates to either target and simulates.
  Through a log or logit link, the output changes materially. For `gaussian` and `exgaussian`,
  only the intercept of the output changes, and no power or precision result for a coefficient
  does. `scale_intercept = TRUE` restores the old behaviour.

* Each point of a sweep draws its own replicates. Every grid point reused the specification's
  seed, so neighbouring points shared their random numbers, while `solve_curve()` treats them as
  independent. Over 1,000 emulated curves a setting, its nominal 95% intervals covered the true
  sample size in 88% of solves on a tight grid, and the true effect in 70% on an effect-size
  axis. Grid point i of `sweep_spec()` now takes `replicate_seeds(seed, length(values))[i]`,
  which brings the two to 96.5% and 95.7%. Every number a sweep produces changes, in
  `power_curve_mixed()`, `precision_curve()`, the power curve of the no-code app, installed or in
  the browser, and the HPC precision-array job, whose tasks a local `precision_curve()` over the
  same grid now reproduces. `seeds = "shared"` (`.seeds` in `sweep_spec()`) keeps the old
  behaviour. The coverage figures quoted for `solve_curve()` come from a new calibration,
  `tools/calibration/solve_curve_coverage.R`, where the earlier 35 in 36 rested on 12 curves
  whose points shared their seeds. Rerun with independent points, those 12 curves no longer
  favour the probit link over the logit, and the documentation now says so.

## New features

* `calibrate_response()` gains `denominator = "random"`, which calibrates the random components
  plus the residual. The coefficient of a contrast coded `c(-0.5, 0.5)` is then Cohen's d in the
  sense of Westfall, Kenny and Judd (2014). The documentation no longer calls the total, which
  includes the fixed effects, the right denominator for an effect size, and
  `?response_variance` shows how to turn a published d into a coefficient. In
  `reading_time_continuous`, whose fixed effects carry 15% of the total, `cond` calibrates to
  0.127 against the random components and the residual, and to 0.117 against the total.

* Every rate carries its Monte Carlo error. `power_design()` gains `power_mcse`, `power_lo` and
  `power_hi`, the standard error and 95% Wilson interval of power, which `power_mixed()` already
  reported, and it still returns a plain list. Both functions gain `type_m_mcse`, the standard
  error of Type M, and a Wilson interval for Type S, `type_s_lo` and `type_s_hi`, since the
  standard error of Type S is zero whenever no sign is wrong. `power_curve_mixed()` and
  `sweep_spec()` return them as columns. Type S and Type M rest on the significant replicates
  alone, which are often a handful, yet they were reported without any uncertainty. Printing a
  `power_mixed()` result shows them in a second table, and the no-code app, installed or in the
  browser, shows each rate with its error. The Python twin computes these values with the same
  operations in the same order, and both test suites pin them bit for bit.

* `power_mixed()` and `precision_design()` simulate a preregistered simplification rule. An
  analysis plan that drops the random correlations, then the slopes, whenever a fit is singular
  analyses each replicate with a model that depends on its data. Its power and Type I error rate
  are therefore not those of the maximal model, which pilotr fitted to every replicate. `formula`
  takes a ladder of models, and `fallback`, `"singular"` or `"warning"`, says when a replicate
  moves on to the next rung. Each replicate keeps the first fit the rule accepts, or else the last
  rung's fit. `model_formula()` gains `structure = c("maximal", "zcp", "intercepts")`. A ladder
  written as those names, such as `c("maximal", "zcp", "intercepts")`, is fitted with the model
  the family needs, and a list of formulas is fitted as written. `power_mixed()` returns and
  prints the share of the replicates analysed at each rung as `rung_shares`, and
  `precision_design()`, `precision_curve()` and the sweeps return it as `rung_<name>` columns.
  `precision_curve()` gains `fallback` as well. The power vignette runs the rule on its priming
  design, null condition included, since a rule that can end on random intercepts alone can
  reject a true null hypothesis too often (Barr et al., 2013; Matuschek et al., 2017). A single
  formula, the default included, gives the results it gave before.

* `spec_from_model()` reports which fixed effects it took from the pilot, and gains
  `effects = "safeguard"`, which replaces each effect by the bound of its 60% interval nearer zero
  (Perugini, Gallucci and Costantini, 2014), and `sesoi`, which sets named effects to a smallest
  effect size of interest. A pilot's point estimates are noisy and, when a pilot is followed up
  because it looked promising, biased upwards, so a study powered for them tends to be
  underpowered (Albers and Lakens, 2018). An interval that includes zero leaves even the sign of
  the effect open, so the safeguard then stops with an error that points to `sesoi`. The result
  records the choice in its `effects` attribute. By default, the coefficients are the pilot's
  estimates, as before. `?spec_from_model` explains when each choice fits, and that a named
  `focal` in `power_mixed()` leaves the simulated effect as it is.

* `counterbalance` on a within factor gives each subject one condition per item, rotating the
  lists across subjects and items, as in a Latin-square design. Planning such a study with a
  fully crossed specification doubled the observations per participant and overstated power.
  Subject s meets item t in condition ((s − 1) + (t − 1)) mod C, where C counts the conditions,
  or the combinations of levels when several factors are counterbalanced, and nothing is drawn.
  The field belongs to specification version 0.4, needs an item unit and is refused with
  `per_subject`. A warning says when the lists are unequal, because the subjects or the items are
  not a multiple of C. `model_formula()` and `model_data()` are unchanged, since with at least C
  subjects and C items every subject meets every condition across items and every item every
  condition across subjects. `build_spec()` gains `counterbalance = TRUE`, and the no-code app,
  installed or in the browser, offers it as "Counterbalanced lists" for a within design with
  items. `spec_from_model()` sets the field when it reads a counterbalanced pilot, as the next
  entry describes. The two-list encoding of the specification suits two conditions only, and it
  stays documented for releases that read no version above 0.3.

* The simulated data of a counterbalanced specification name each subject's list in a column
  `list`, placed after `subject` and `item`. Subject s follows list ((s − 1) mod C) + 1, and
  nothing is drawn for it, so a specification without a counterbalanced factor simulates as
  before. Such a design keeps the name `list` for that column, and `validate_spec()` refuses a
  factor, predictor, grouping factor, contrast column or response of that name in it.
  `spec_from_model()` recognises a counterbalanced pilot and returns a specification-0.4 design
  with `counterbalance` set, and says so in a message. `counterbalance = FALSE` keeps the earlier
  output, a crossed specification that releases reading 0.3 accept, with a warning. A pilot whose
  subjects each saw only some of the items keeps `per_subject` and stays crossed, with the
  warning, since the field needs every subject to see every item. A pilot whose model has a term
  named `list`, as when the list is fitted beside the condition, stays crossed too, with a
  warning that asks for that term to be renamed.

* `power_design()`, `power_mixed()` and `power_curve_mixed()` gain `alternative`, which plans
  the one-sided test of a directional prediction with `"greater"` or `"less"`. Every earlier
  result was two-sided, and the default, `"two.sided"`, leaves every number as it was. The
  one-sided p-value comes from the statistic behind the two-sided one, a t with its degrees of
  freedom from `t.test()`, `lm()` or lmerTest, or a Wald z. Type S and Type M are taken over the
  replicates significant in the tested direction. `power_design()` records the test in
  a new element, `alternative`, and the `fitter` of `power_mixed()` names the direction. On the
  same replicates, a one-sided test at alpha rejects exactly where the two-sided test at
  2 alpha rejects with the estimate in the tested direction. A two-sided result at 2 alpha from
  an earlier version therefore gives the one-sided power at alpha, in the direction of the true
  effect, as `power * (1 - type_s)`.

* `focal` in `power_mixed()` and `power_curve_mixed()` accepts named sets of coefficients, each
  tested jointly, so the power of a multi-degree-of-freedom effect, such as a three-level factor
  or its interaction, can be planned. `list(interaction = c("A1_B", "A2_B"))` tests the two
  coefficients of a 3 × 2 interaction together. Each coefficient could only be tested by itself,
  and the power of the joint test differs from that of either coefficient. A linear mixed model
  has Satterthwaite's F from `lmerTest::contest()`, a `glmer()` or `glm()` fit a Wald chi-square
  and an `lm()` fit the F test of the model with and without the set, so a set of one
  coefficient reproduces the single-coefficient test. Each set's power carries its Monte Carlo
  error and Wilson interval, and `fitter` names the joint test. A joint test has neither a single
  estimate nor a direction, so Type S and Type M are `NA` and a one-sided `alternative` is
  refused. The result records the sets, its printout lists their coefficients and a sweep adds a
  `coefficients` column. `precision_design()` and `precision_curve()` refuse sets, since a set has
  no single interval. The Python twin's `power_mixed()` has no joint test, and its documentation
  says so.

* `power_report()` writes a sample-size justification in Markdown from a power result and,
  optionally, a solved curve. The documentation called the sample size from `target_n()` the
  number a preregistration quotes, but the justification behind it had to be assembled by hand.
  The report gives the items Lakens (2022, Table 4) asks for. It states the primary analysis with
  its fitter, test, alpha and sidedness, and the desired power with the solved sample size, its
  interval, transform and heterogeneity factor. Each focal effect follows with the justification
  given for its size, and `null_result` gives the power when every focal effect is zero. The
  simulation is accounted for as Morris, White and Crowther (2019) recommend, with the
  replicates, the fit counts, the share analysed at each rung of a simplification rule and the
  Monte Carlo standard error and Wilson interval of every rate, computed in the report from the
  counts. The report ends with the versions of pilotr, R, lme4 and lmerTest, the specification's
  seed and an MD5 checksum of the specification file's bytes. For a `precision_design()` result, each focal effect has its
  mean interval width and the probabilities of the two decisions about the region of practical
  equivalence, which `rope` states. The Python twin's `power_report()` writes the same text, byte
  for byte, from the same counts, which `tools/parity/report_cross.py` checks in CI. The power
  vignette writes one for its priming design.

## Fixes

* Printing a `power_mixed()` result applies `digits` to the interval bounds as well as to the
  rates. The bounds were written with three decimals whatever `digits` said, and the
  documentation called `digits` significant digits, where it counts decimal places.

* The no-code app, installed or in the browser, sends a design its power tab cannot analyse to
  the analysis that fits it. It sent every such design to the `power_mixed()` of both packages.
  The Python one is now offered only for one within factor crossed with items, in a family whose
  coefficients are on the scale it analyses, and ordinal and Beta designs are sent to
  `generate_design_analysis()` and `brms_bridge()`.

* The no-code app, installed or in the browser, keeps the session open after an input its power
  tab cannot use, and says what is wrong. An N of 2, a cleared N subjects box or a pasted
  three-level design ended the session, and the design was lost with it. The Power tab now asks
  for at least 4 subjects and reports any other error as text. In the installed app, a pasted
  specification that `load_spec()` refuses made the power buttons do nothing, and its refusal now
  appears there.

* In the no-code app, installed or in the browser, the Data tab, the summary, the plot and the
  data download show the last simulation, labelled with its own seed, and a note above the tabs
  says when the controls have changed since. The summary and the plot failed after a change of
  family or response name and lost the grouping after a change of factor name. The data download
  was simulated afresh from the controls, so it could differ from the table on screen, and the
  installed app's Data tab printed the Seed box, which after a change or under a pasted
  specification was not the seed of the data. The download is now named after the simulated
  specification and its seed. The plot of a pasted design whose factor has more than two levels
  used to stop with "Insufficient values in manual scale" and now draws.

* The no-code app's power curve, installed or in the browser, extends its grid until it brackets
  the 0.80 target. The grid stopped at twice the current N, which for the default design is the
  target itself (127.5 subjects by `power.t.test()`), so about half of all runs were refused and
  those that solved were biased low, at 122 subjects on average in an emulation of 3,000 runs. A
  point half as large again is now added above the grid while the top point's power falls short
  of 0.80 by more than twice its Monte Carlo error, and a point two thirds the size is added below
  it while the bottom point's power exceeds 0.80 by as much. The grid stops at 8 times N in the
  installed app and at 2,000 subjects in the browser, and holds 10 points at most. Each point keeps
  the seed it was first given. At the six seeds the app's tests use, the default design now solves
  every time, to a mean of 129 subjects, where three of the six were refused. For a design without
  random effects or continuous predictors, the Power tab also gives the sample size
  `power.t.test()` computes, as a check on the curve's.

* `run_app()` restores the caller's future plan when the app stops, sets a plan only over a
  sequential one, with two workers in place of one for every core, and honours `async = FALSE`,
  which the app ignored. It used to end by setting a sequential plan, discarding any plan set
  before. With `async`, the power curve now runs in the background process as well as the power
  analysis. The documentation of `run_app()`, the README and the article on the no-code app say
  what this does: it frees the app's R process for other sessions, while the session that started
  a run waits for its result. They used to say that the app does not block and that its
  simulations use all of your local cores.

* `n_warning` no longer counts boundary-singular fits, which `n_singular` already counts. lme4
  records its singular-fit notice among the optimiser's messages, so every singular fit was also
  counted as a fit with warnings, and a run could not show whether any fit had failed to
  converge. `n_converged` and `n_singular` are unchanged.

* Every specification field is read by its exact name. R looked fields up with `$`, which
  matches a name partially, so a design whose random entries were `subject_site` and no
  `subject` acquired by-subject random effects it never declared, and its data differed from the
  Python twin's in every row. With `strict = FALSE`, a misspelt field such as `rounding` acted as
  the field it resembled, although the warning said it would be ignored.

* `load_spec()` and `validate_spec()` refuse a specification that repeats a key within one
  object. JSON leaves a repeated key undefined and the two readers disagree: `{"grp": 0.3,
  "grp": 0.5}` applied 0.3 twice in R and 0.5 once in Python. The refusal does not depend on
  `validate`.

* A seed must lie within ±(2^53 − 1), the range in which both twins read a JSON integer exactly.
  Larger seeds were read as different numbers by the two readers, which then drew different data.

* `build_spec()` and `spec_from_model()` keep a seed of 2^31 or more, which they used to turn
  into `NA`, and leave a non-whole seed for `validate_spec()` to refuse, where they used to
  truncate it. Seeds in the integer range are stored as before. The app accepts a larger seed
  for the same reason and prints it in full.

* A factor level, a `vary_within` entry or a `correlated` flag written as `null` is refused, as
  the Python twin refuses it. R labelled those rows `NA` and read the flag as `false`.

* A malformed `response`, `predictors`, `random` or `spec_version` value is reported in the
  validator's own words. Some shapes stopped with a base R error, and a `spec_version` given as
  an array of two versions was read as the first.

* `load_spec()` reads a file that starts with a byte-order mark without warning about it. The
  specification format allows one, since some Windows editors add it.

* `design_conditions()` changes only the effects it names. Through `sweep_spec()`, its
  conditions replaced the whole coefficient list, so every coefficient they did not name was
  dropped from data generation, from `model_formula()` and from the default focal set, contrary
  to the documentation. On `reading_time_continuous`, a sweep of `cond` simulated and analysed a
  one-coefficient model. Conditions are now merged into the specification's own coefficients, in
  its order, and a condition naming a coefficient the specification lacks is refused. They carry
  the class `pilotr_conditions`, which subsetting them with `[` and joining them with `c()` keep.
  Any other list still replaces the addressed field wholesale. With 0.3.x, pass
  `.base = spec$fixed$coefficients` to `design_conditions()` for the same result.

* The `sweep_spec()` example no longer runs the null condition twice.

* `pilotr::run_app()` works in a session that has not attached pilotr. The app looked for the
  package's functions among attached packages and then for source files, and in an installed
  package found neither, so the documented launcher stopped with "no pilotr engine sources". The
  app now takes the functions from the package's namespace, without attaching the package.

* `validate_spec()` refuses specifications whose names collide. Each of these validated and then
  moved an effect, rescaled a variance or overwrote a column without a word: two columns with one
  name (a factor called `subject`, a response named like a factor), a contrast column defined by
  two factors or named like a predictor, a unit, a grouping factor, another factor or the
  response, an interaction whose analysis column `a_b` already exists, a level listed twice, a
  correlation pairing a term with itself or giving one pair twice, a factor both between and
  within a unit, and a `random.item` entry in a design without items, which was dropped. Each
  refusal reads the same in the Python twin. The rules are set out under Names in the
  specification.

* `validate_spec()` refuses a factor, predictor or extra grouping-factor name that is empty or
  made only of spaces, and a response name made only of spaces (an empty one was already
  refused). Such a name validated and then became a column with no visible name. An empty factor
  name stopped `simulate_design()` with base R's "replacement has length zero", which does not
  say which field is at fault, and a blank `random` key was reported as
  `random. must be an object`.

* Two refusals read as the Python twin's do. A `correlated` flag that is not a boolean must be
  `true or false`, as the JSON file spells it, and `response.thresholds` must be "a number or a
  non-empty numeric array".

* `model_data()` no longer fills a factor's later contrast columns with `NA` when the first
  contrast column carries the factor's own name. With three levels, every `power_mixed()`
  replicate then failed to fit.

* `power_design()` reports a true effect of 0, with Type S and Type M undefined (`NaN`), when the
  specification has no coefficient for the factor's contrast, as with `"coefficients": {}`. It
  stopped with "missing value where TRUE/FALSE needed".

* `precision_design()` and `precision_curve()` check `rope`. A vector was recycled across
  replicates, and a negative value was used as given. On the design of the `precision_design()`
  example, the range form `c(-0.02, 0.02)` reported a `p_meaningful` of 0.200 over 30
  replicates, where `0.02` gave 0.067. `rope` takes a half-width or an increasing pair
  `c(lower, upper)`, which also allows an asymmetric region, and anything else is refused.
  Results for a non-negative half-width are unchanged. The HPC precision-array job reads its
  `ROPE` variable in either form, as `0.05` or `-0.02,0.05`, and checks it in the same way.

* A Poisson mean above about 708 is refused. The sampler's starting value `exp(-mean)` is
  subnormal there, and from a mean of about 738 the counts were biased or intermittently refused,
  while the specification promised service up to 746. Over 40,000 draws, a mean of 744 drew counts
  with a mean of 733.5 and a variance of 407, and a mean of 742 was refused for 4% of them. The
  refusal reads as the Python twin's and asks for a mean of at most 708.

* `spec_json()` writes each number so that 'jsonlite' and 'Python' read back exactly the double
  R holds. Version 0.3.1 checked each candidate width with `as.numeric()`. On builds with
  long-double arithmetic, among them R on x86 Windows and Linux, `as.numeric()` rounds a decimal
  twice, so the check accepted some 15- and 16-digit forms that a correctly rounded reader puts
  one unit in the last place away, about one number in 5,000 to 7,000. A specification saved and
  reloaded then simulated different unrounded data, in R and in 'Python' alike. Candidates are
  now checked with 'jsonlite', the reader `load_spec()` uses.

* `generate_r_script()` and `generate_design_analysis()` embed the specification as the JSON
  text `spec_json()` writes, and the script reads it with `load_spec()`, so 'jsonlite' recovers
  every number exactly on every build of R. The specification used to be written as R code. A
  build without long-double arithmetic, such as R on Apple silicon, can parse a 17-digit decimal
  in R code one unit in the last place away, and likewise a 16-digit one whose digits exceed
  2^53. A script carrying a fitted coefficient could then simulate different data there. In a
  design analysis, the decision rule, the gate and the true values stay plain decimals wherever
  every build reads them exactly, as it does any value typed with up to four decimal places. Any
  other value is read through 'jsonlite'.

* Beta responses stay inside (0, 1). Near a ceiling, `x / (x + y)` rounds to exactly 1 in double
  precision (161 of 4,000 rows at an intercept of 3 and `phi = 2`), and rounding sent low
  proportions to exactly 0 (17% of the rows at three decimals in one realistic design). brms's
  `Beta()` refuses either value. An exact 1 or 0 now becomes the nearest double inside the
  interval, a rounded value at 0 or 1 becomes the nearest value inside it at that many decimals,
  and `round = 0` is refused for Beta. No other value changes, so the shipped examples simulate
  the same data as before.

* Extreme linear predictors are refused with the same message in both twins. R returned `Inf` for
  a lognormal or shifted lognormal response whose `exp(eta + sigma * z)` overflowed, and 0, 1 or
  `NaN` for a Beta response whose shape parameter was zero or whose two Gamma variates both
  underflowed. Python raised a bare `OverflowError` or `ZeroDivisionError` in the same cases.

* `power_curve_mixed()` and `precision_curve()` start one worker pool for the whole sweep, as their
  documentation said. Since 0.3.0 every grid point started and stopped its own, at a cost of 20 to
  50 seconds a point on a loaded machine. `power_mixed()`, `precision_design()` and
  `power_design()` gain `cl`, for a cluster the caller has made, which they run on and leave
  running. Passed through `sweep_spec()`, as in `sweep_spec(spec, path, values, power_mixed,
  cl = cl)`, one cluster serves every point of an effect-size sweep as well.

* The no-code app's Verify button runs the script that `generate_r_script()` exports, in a clean
  R session, and reports reproduction only when the script's data are identical, column by
  column, to the data the current specification produces. Otherwise, it names the first column
  and row that differ. It used to simulate from the JSON in that session without running the
  script, and compared the sums of the response column with a tolerance, so it certified a script
  it had never run. A specification the engine refuses is now reported beside the button, where a
  cleared N subjects box used to end the session.

* Exported scripts give the CRAN installation line, naming lme4 and lmerTest only where the
  suggested analysis needs them, and suggest the analysis that fits the design:
  `power_design()` for a two-group Gaussian design with one row per subject,
  `generate_design_analysis()` for ordinal and Beta outcomes, and `power_mixed()` for any other
  design. Every script used to suggest `power_mixed()`, and its installation line said "once
  available", although pilotr is on CRAN.

* `generate_r_script()` reads a specification given as the path of its file. It used to embed
  the path itself, so the script ran only while the file stayed at that path.

* The tests need testthat 3.1.5, which added `expect_no_error()` and `expect_no_warning()`, so
  DESCRIPTION now asks for it in place of 3.0.0. The continuous-integration job meant to test the
  declared floors had never installed them: on R 4.6, testthat 3.0.0 does not compile and
  ggplot2 3.4.0 does not load, and the job went on to pass on the current releases. It now runs
  on an older R, reads the floors from DESCRIPTION and fails unless each one is installed.

* The browser build of the no-code app gives the command that installs the Python package as
  the bare `pip install "pilotr[mixed]"`, with its explanation on the line before. The command
  ended with a `#` comment, which an interactive `zsh`, the default shell on macOS, passes on to
  pip, so a pasted line stopped with "Invalid requirement: '#'". The R installation line keeps
  its trailing comment, which R ignores.

## Deprecated

* A within factor whose `vary_within` omits a unit of the design now warns, in strict and lenient
  validation alike. pilotr has always crossed a within factor with every unit, so the list
  changed nothing. From spec version 0.4, such a list is refused.

## Documentation

* `citation("pilotr")` names the R package version alone. It read "R and Python package
  version 0.3.1" while the Python package stood at 0.3.0, and the two packages are released
  separately. It also gives the author's ORCID iD now.

* `citation("pilotr")` and the website give the DOI that CRAN assigned to the package,
  <https://doi.org/10.32614/CRAN.package.pilotr>, in place of the Zenodo concept DOI. It names
  the R package across all its versions, and it is the DOI that the CRAN page and the citation
  R generates for a CRAN package already gave, so the R package no longer circulates under two.
  `CITATION.cff` describes both packages, so it keeps the Zenodo DOI, which archives every
  GitHub release of both alike, and now lists each package's identifiers, marked with the
  package they name. The repository's `README.md` gives a reference for each package.
  The CRAN DOI and the Zenodo DOI are shown as full links on the website and in the `README.md`
  files. The R package's `README.md` gains a badge for the DOI CRAN assigned, for GitHub and
  CRAN, which the home page leaves out of its sidebar, where the DOI is already given under
  Citation.

* The 0.3.0 and 0.3.1 entries below no longer give release dates. Both gave dates on which
  neither version reached CRAN, and the site already prints the CRAN date beside each entry.

* The tree CRAN published as pilotr 0.3.0 is tagged `cran-0.3.0` on GitHub, and 0.3.1 is
  tagged `v0.3.1` and archived on Zenodo. The `v0.3.0` tag and its Zenodo record predate two
  changes that CRAN's 0.3.0 carried, namely that `brms_bridge()` returns a visible
  `pilotr_bridge` object and that this object has a print method.

* `sweep_spec()` and `solve_curve()` say that a sweep over lists, `design_conditions()`
  included, records a grid index in its leading column, and that `solve_curve()` solves for that
  index unless `x` names another column. Both pages show the sweep to use for a single effect,
  `sweep_spec(spec, "fixed$coefficients$cond", values, power_mixed)`. It keeps the other
  coefficients, and `solve_curve()` returns a minimum detectable effect from its curve.

* The specification states how units are allocated: the block rule for between factors and
  extra grouping factors, the nesting of a between factor within clusters that follows from it,
  and the full crossing of within factors. It gives worked encodings for a 2 × 2 between design,
  a two-list counterbalanced design and randomisation within clusters, each held bit-identical
  across the twins by the parity harness. The power vignette says that its priming example is
  fully crossed, and that simulated data are complete, so expected attrition is allowed for by
  recruiting N / (1 − p).

* The precision vignette says that its priming example is fully crossed, and that a
  counterbalanced study has half the observations and is planned with `counterbalance = TRUE` in
  `build_spec()` or with the two-list encoding. `?spec_from_model` explains the counterbalancing
  warning, the field and that encoding. It also shows a `prep` function for `power_mixed()` and
  `precision_design()` that drops the rows a counterbalanced study would not observe and runs in
  the worker processes of `workers > 1`.

* The cross-language guarantee is stated as it holds, in the specification, the README files,
  `DESCRIPTION`, the vignettes, the apps' guides and the landing page. R and Python produce
  bit-identical data wherever they share a maths library, as on Linux. On Windows they do not,
  and a value that passes through `exp()`, `log()` or `pow()` can differ in its last bits.
  Measured there, about 1.4e-5 of normal draws differed, since the tail branch of AS 241 calls
  `log()`, so Gaussian data are not exempt. Lognormal and ex-Gaussian responses differed by one
  unit in the last place and unrounded Beta responses by up to 35, and a discrete outcome can in
  principle move by a whole category. The guarantee used to cover the Gaussian family everywhere
  and every family with `round` set, although `round` is refused for the discrete families.
  Rounding has removed every continuous difference observed, but it guarantees nothing. One Beta
  case can differ on Linux too. From R 4.6.0, R's `^` computes a cube or a fourth power by
  repeated multiplication, where Python calls `pow()`. A Beta response whose Gamma shape is
  exactly 1/3 or 1/4 can therefore differ in its last bits on any platform, and R's own data for
  such a design changed with R 4.6.0. Within one language, the data are identical on one
  platform, and across operating systems on the same terms as across languages. The vignettes
  show determinism with `identical()`, since `all.equal()` allows a relative difference of 1.5e-8.

* The specification documents every sampling rule a third implementation needs, in a section of
  pseudocode: the partial Fisher-Yates shuffle, the block rules, uniform predictors and
  reliability, the Cholesky factorisation, the two-branch inverse logit, each family's response,
  the Poisson walk, the ordinal rule, the Marsaglia-Tsang Gamma sampler and the Beta ratio. It
  corrects the Gamma sampler's draw consumption. A proposal draws a uniform only when its normal
  gives a positive `v`, and a boosted variate draws its uniform after the Gamma(shape + 1)
  variate. The schema refuses an unknown top-level field and a response parameter its family
  does not use, requires each family's parameters, and accepts the single-value forms of
  `vary_within` and `thresholds`, as strict validation does. `?build_spec` says why a Beta
  response is left unrounded: at four decimals, rounding put 9.8% of the rows of one design on
  0.0001 or 0.9999, and such rows dominate a Beta likelihood.

* `?replicate_seeds` gives the right reason why consecutive seeds correlate. Seeding is affine in
  the seed, so the first draws of consecutive seeds step by a constant mod 1 and correlate 0.95,
  and a longer warm-up changes only the step (Matsumoto et al., 2007). It said that consecutive
  replicates began a few steps apart in one sequence.

* The power and precision vignettes show results produced by this version, with the denominator
  of power stated correctly (`n_returned`), Wilson intervals taken from the returned columns and
  the model that `power_mixed()` fits. Their cached results predated 0.3.0. The precision
  vignette no longer says that every replicate converged, and it computes each figure it quotes
  from the cached sweep. `tools/vignette-caches.R` in the source repository regenerates the
  caches, and a test fails when the seeds, the generative engine or a cached chunk change and the
  caches do not.

* `?power_curve_mixed`, `?solve_curve` and the power vignette say that when an effect varies from
  item to item, power at a fixed number of items levels off below 1, however many subjects are
  added (Westfall, Kenny and Judd, 2014; Judd, Westfall and Kenny, 2017). A curve over subjects
  that levels off near the target can then be solved far from the true crossing, or as though it
  crossed when it never does. The vignette shows how to compare designs that trade subjects for
  items, by sweeping the number of items at a few subject counts and solving each curve.

* The HPC job runs on any SLURM cluster and calls `precision_design()`, so its output has the
  same columns and fit accounting as `precision_curve()`. The job charged the author's account,
  wrote under the author's project paths and ran a copy of its runner kept outside the
  repository, with a bootstrap script the repository did not hold. Its runner repeated the
  replicate loop and counted every fit that returned an estimate as converged, the definition
  0.3.0 replaced. The wrapper now marks the account and the project directory as the two values
  to fill in, and runs the runner saved beside it. The README gives the commands that install
  pilotr, lme4 and lmerTest into the project's library, shows how the task files combine into a
  curve for `solve_curve()`, and limits the claim of exact reproduction to the simulated data,
  since the model fits depend on the lme4 release and the platform.

* The power vignette's table of fitters, the cross-language vignette and `README.md` describe the
  Python twin's `power_mixed()` as it now stands. It fits the model the specification implies,
  with each random effect as an independent variance component, from ports of `model_formula()`
  and `model_data()`. `tools/parity/formula_cross.py` holds the ports to this package's functions
  on every shipped example and parity case, the formulas character for character and the contrast
  and product columns bit for bit.

* Both apps' guides say that crossed mixed-effects power runs in R through lme4, the reference
  backend, and in Python through statsmodels. They named lme4 for both. Before a run, the
  installed app's power tab said that simulations were unlimited once you install the package,
  although the app is part of the installed package. It now gives the app's cap, which
  `PILOTR_MAX_SIMS` sets, and says to call `power_design()` directly for more. The vignette on
  the no-code app names the tabs of both apps. It named the browser build's alone.

* The website's R pages link favicons that exist. pkgdown 2.1.2 changed the set of favicon files
  its pages link, and three of them, `favicon-96x96.png`, `favicon.svg` and `site.webmanifest`,
  returned 404 on every page. `tools/favicons/build_favicons.R` in the source repository builds
  the set offline. In a browser tab, the browser build of the app is titled "pilotr: no-code
  app", not "Shiny App", and shows the pilotr icon. The theme toggle on the site's landing and
  not-found pages works in a browser that blocks site storage, where its script used to stop
  before the toggle was wired up. The contributing guide names the checks in `tools/parity/` that
  CI runs, in place of a comparison of CSV files to a tolerance of 1e-6. The note on advanced
  designs in the repository describes specification 0.3.

## References

Albers, C. and Lakens, D. (2018). When power analyses based on pilot data are biased: Inaccurate
effect size estimators and follow-up bias. *Journal of Experimental Social Psychology, 74*,
187-195. <doi:10.1016/j.jesp.2017.09.004>

Barr, D. J., Levy, R., Scheepers, C. and Tily, H. J. (2013). Random effects structure for
confirmatory hypothesis testing: Keep it maximal. *Journal of Memory and Language, 68*(3),
255-278. <doi:10.1016/j.jml.2012.11.001>

Bolker, B. M., Brooks, M. E., Clark, C. J., Geange, S. W., Poulsen, J. R., Stevens, M. H. H. and
White, J.-S. S. (2009). Generalized linear mixed models: A practical guide for ecology and
evolution. *Trends in Ecology & Evolution, 24*(3), 127-135. <doi:10.1016/j.tree.2008.10.008>

Judd, C. M., Westfall, J. and Kenny, D. A. (2017). Experiments with more than one random factor:
Designs, analytic models, and statistical power. *Annual Review of Psychology, 68*, 601-625.
<doi:10.1146/annurev-psych-122414-033702>

Lakens, D. (2022). Sample size justification. *Collabra: Psychology, 8*(1), 33267.
<doi:10.1525/collabra.33267>

Li, P. and Redden, D. T. (2015). Comparing denominator degrees of freedom approximations for the
generalized linear mixed model in analyzing binary outcome in small sample cluster-randomized
trials. *BMC Medical Research Methodology, 15*, 38. <doi:10.1186/s12874-015-0026-x>

Matsumoto, M., Wada, I., Kuramoto, A. and Ashihara, H. (2007). Common defects in initialization
of pseudorandom number generators. *ACM Transactions on Modeling and Computer Simulation,
17*(4), 15. <doi:10.1145/1276927.1276928>

Matuschek, H., Kliegl, R., Vasishth, S., Baayen, H. and Bates, D. (2017). Balancing Type I error
and power in linear mixed models. *Journal of Memory and Language, 94*, 305-315.
<doi:10.1016/j.jml.2017.01.001>

Morris, T. P., White, I. R. and Crowther, M. J. (2019). Using simulation studies to evaluate
statistical methods. *Statistics in Medicine, 38*(11), 2074-2102. <doi:10.1002/sim.8086>

Perugini, M., Gallucci, M. and Costantini, G. (2014). Safeguard power as a protection against
imprecise power estimates. *Perspectives on Psychological Science, 9*(3), 319-332.
<doi:10.1177/1745691614528519>

Westfall, J., Kenny, D. A. and Judd, C. M. (2014). Statistical power and optimal design in
experiments in which samples of participants respond to samples of stimuli. *Journal of
Experimental Psychology: General, 143*(5), 2020-2045. <doi:10.1037/xge0000014>

# pilotr 0.3.1

## Fixes

* `spec_json()` no longer corrupts a coefficient on a build of R without long-double
  arithmetic. Numbers were written by asking `jsonlite` for 17 significant digits and then
  shortening each one in the text, which meant reading the number back with `as.numeric()`.
  That read-back is not correct everywhere. Without long doubles `as.numeric()` accumulates
  the mantissa in a double before applying the decimal exponent, and 17 significant digits
  overflow the 53-bit mantissa on the way, so `"0.33333333333333331"` came back one unit in
  the last place below `1/3`. The shortened form then recorded that wrong value, and the
  saved specification no longer held the coefficient the user had set. Each number is now
  formatted from the double itself, before `jsonlite` sees the document, so nothing is ever
  read back from the text it was just written to. A coefficient of `1/3` is written at 16
  significant digits rather than 17, which is both the shorter form and the one such a build
  reads back exactly. Found by CRAN's no-long-double check.

* `generate_r_script()` could write `Inf` in place of a finite coefficient on the same builds.
  It embeds each number through the helper that looks for the shortest decimal that reads back
  as the same double, and that helper returned nothing at all when none of the three widths it
  tried passed the check, a case the caller took to mean the value was infinite. Without long
  doubles the 17-digit width fails that check for the very values that need it. The helper now
  falls back to 17 significant digits, which is the widest form there is and so the one to use
  when nothing shorter will do, and returns nothing only for a value that really is not finite.

# pilotr 0.3.0

## Read this first

Two changes alter numbers that earlier versions produced. Both are deliberate corrections, and
both are unavoidable if the results are to mean what the documentation says they mean. Install
0.2.1 to reproduce output from 0.2.1.

* Replicate seeds changed. The power and precision loops previously seeded replicate `i` with
  `seed + (i - 1)`. Consecutive seeds are not independent streams in this generator, and the
  measured consequence was severe. The first draw of replicate `i` correlated 0.95 with the first
  draw of replicate `i + 1`. Seeds are now drawn from the shared generator, which brings that
  correlation to -0.02 and moves a Ljung-Box test over the replicate means from p below 0.0001 to p
  of 0.94. Every number the replicate loops produce is therefore different, and the first replicate
  no longer uses the specification's own seed. `replicate_seeds()` is exported so that a
  hand-written loop or a cluster array task can use the same rule.

* An interaction random slope now reaches the data. A slope keyed on an interaction, such as
  `"z_cosine:z_ISI"`, was accepted, sized into the covariance, drawn from the random stream and
  then multiplied by a lookup that returned nothing, so it was silently discarded. Meanwhile
  `model_formula()` and `brms_bridge()` both emitted that slope, so the model fitted was richer than
  the process that generated the data, which inflates power in the direction Barr et al. (2013)
  describe. Any design that used one was anti-conservative and its data have changed.

## Cross-language reproducibility

The R and 'Python' implementations agreed on far less than the documentation claimed. Four separate
defects made them diverge, none of them visible in the shipped examples, because six of the eight
set `response.round`, which quantises exactly these differences away. With rounding removed, three
of the eight examples differed, in up to 9.24% of cells.

* The accumulation order of the linear predictor is now fixed by the specification. R folded the
  terms one at a time while 'Python' summed them and added the total. Floating-point addition is not
  associative, so the two disagreed on 63.8% of rows in a matched reproduction over 200,000 rows.
* Neither language's built-in summation is used any more. Base R's `sum()` accumulates in 80-bit
  long double, and CPython's `sum()` has applied Neumaier compensation to floats since 3.12. Both
  are more accurate than a plain double fold, in different ways, so an inner product of three terms
  or more could land on different doubles. Every inner product is now an explicit fold.
* Integer powers are written as repeated multiplication. R special-cases small integer
  exponents while 'Python' calls the library `pow()`. Over 200,000 draws in the Gamma sampler's
  range the two disagreed on a third of inputs by up to 6 ulp. Because that value decides a
  rejection step, it also changed how many draws were consumed. This alone accounted for every
  difference in the `beta_proportion` example.
* The random-effect covariance is bracketed identically. R computed `(sd_i · sd_j) · r_ij` and
  'Python' `(sd_i · r_ij) · sd_j`. Multiplication is commutative but not associative, and the
  difference propagated through the Cholesky factor into every random effect drawn.

All eight shipped examples are now bit-identical between the two languages as shipped, as are five
adversarial specifications added to exercise what the shipped examples do not. With `response.round`
removed, six of the eight remain bit-identical, and the two whose family applies `exp()` to the
linear predictor, `crossed_mixed_rt` and `reading_time_continuous`, differ by at most one unit in the
last place, on 0.09% and 0.05% of cells. That residue is the libm limit described immediately
below, and `tools/parity/tolerance.json` records which cases carry an allowance for it, and why.

* The scope of the guarantee is now stated honestly in `spec/SPEC.md`. It is exact for
  `gaussian`, for any design applying no transcendental function to the linear predictor, and for
  any family with `response.round` set. The families that apply `exp()` or `log()` may differ in the
  last unit in the last place, because IEEE-754 does not require correct rounding for those
  functions and the two builds need not share a maths library. Measured rates are given, and the
  attribution is demonstrated: the same design switched to `gaussian`, with an identical seed,
  structure and draw sequence, is bit-identical.

* Three small divergences between the two engines are closed. A non-whole seed truncates
  identically, a whole `spec_version` reads the same however it was written, and a
  non-object unit is now reported as such, where R used to crash with a base error. The
  message text of each is character-for-character identical across the twins.

## Validation and versioning

* `validate_spec()` checks a specification against the schema and against the cross-field rules
  the schema cannot express, and `load_spec()` now calls it by default. A strict draft-07 schema had
  shipped since 0.1 with no code path consulting it.
* Validation exists because several ways of getting a specification wrong produced plausible data
  and no error at all. A mistyped coefficient key resolved to no column and so silently set that
  effect to zero, which generates exactly the data of a null design and reports success. A
  response parameter left over from another family was ignored. Both are now refused.
* `spec_version` negotiates the format. A specification with no such field is a 0.2
  specification. One using a 0.3 feature must declare 0.3, because a 0.2 implementation reads it
  differently and generates different data without complaint. One declaring a version newer than
  the implementation understands is refused outright.
* `simulate_design()` gains `validate`, defaulting to `TRUE`. The replicate loops validate once and
  then skip it, so a sweep pays the cost once.

* `spec_from_model()` gains test coverage: the recovered specification's units,
  `between`/`vary_within` placement, interaction keys read back off product columns,
  and random-effect estimates are checked against the design that generated the pilot
  data, alongside the refusal paths for models the function cannot read.

* The package now declares a minimum R version, `Depends: R (>= 4.0.0)`. The cross-language
  claim assumes the `round()` that arrived in 4.0.0, which is where that floor comes from. A
  new CI job checks the declared Suggests floors, which the matrix (release, devel, oldrel-1)
  sat well above and so never exercised.

* The packaged example specifications are now tested against the repository's own copies.
  `spec/examples/*.json` is canonical and both packages carry a mirror so an installed
  copy can reach it, but nothing enforced the mirror. A load-and-simulate test cannot: a
  stale packaged copy still loads and simulates perfectly well, it simply describes a
  different design from the one the repository documents. The new `test-examples.R`
  compares the bytes, and is twinned with the Python suite's `test_examples.py`. It skips
  when the package is checked in isolation from the repository.

## Fixes

* `n_converged` was not the number of converged fits. The replicate loops wrapped each fit in
  `suppressWarnings()`, so only a hard error was visible and every boundary-singular fit counted as
  converged. In the package's own documented example, 12 subjects by 12 items with a maximal model,
  85 of 100 replicates were boundary-singular while `n_converged` reported 100. The loops now report
  `n_attempted`, `n_returned`, `n_converged`, `n_singular` and `n_warning` separately. Singular and
  warning fits are still used, since their fixed-effect estimates remain interpretable and dropping
  them would bias the result.
* A non-positive-definite random-effect covariance is now an error, naming the grouping factor
  and the column at which the factorisation failed. The failing pivot was previously clamped at zero
  and the result returned in silence, which produced random effects whose standard deviations were
  several times those requested. In one test a requested 0.200 came back as 0.805.
* `spec_json()` lost precision on nearly every number. It wrote through `as.character()`, so a
  coefficient of `1/3` round-tripped to `0.33333333333333298`. Over a sample of 214 doubles, 189
  failed to round-trip. Since the JSON file is the portable artefact, the specification itself was a
  source of divergence. Numbers are now written at the shortest precision that round-trips exactly,
  so `0.3` still reads as `0.3`. The blanket `auto_unbox` is gone too, so a one-element
  `vary_within` or a single ordinal threshold stays an array.
* `generate_r_script()` embedded the specification through `deparse()`, which prints 15
  significant digits and so does not round-trip. It now emits numbers at full precision, which
  matters because the point of the script is bit-for-bit reproduction.
* `brms_bridge()` wrote its model to the console instead of returning it. The code reached
  standard output on every call, whether or not anything was there to read it, and the one internal
  caller, `generate_design_analysis()`, had to divert that output around itself through a sink to a
  temporary file. It now returns a `pilotr_bridge` object, visibly, carrying the same `formula`,
  `family`, `priors` and `code` as before, and the new `print()` method writes the code. A bare call
  at the console still shows the model; an assignment is silent.
* `model_formula()` and `brms_bridge()` emitted correlated random effects unconditionally, along
  with an LKJ prior, while the generative process only correlates them when `correlations` is
  supplied. They now follow the new `correlated` flag and emit a double bar otherwise. A group with
  no slopes keeps a single bar, since `lme4` cannot parse `(1 || g)`.
* `model_data()` did not create the product column an interaction random slope needs, so the
  emitted formula referred to a variable the modelling data lacked. It now covers the union of the
  fixed-coefficient and random-slope keys.
* A focal effect that never appears in any fit now warns. It used to return decision proportions
  of zero, which read as 'this design can decide nothing' when the cause was a name that did not
  match the model.
* A replicate loop in which no fit succeeds now reports why, passing on the fitter's own
  message. An unidentifiable random-effects structure previously produced a silent result of `NA`.
* A correlation naming a random-effect term that does not exist now raises a clear error, where it
  used to fail at a subscript.

* A Poisson mean beyond the sampler's reach was returned as the iteration cap. The
  inverse-CDF walk starts from `exp(-mean)`, which underflows to exactly zero once the
  mean passes about 746 (a poisson intercept of 7 already implies a mean of exp(7), about
  1097), after which every simulated count came back as the cap of 1000000 while reporting
  success. Both engines now refuse such a mean, naming `exp(eta)` and the offending value,
  with message text character-for-character identical across the twins. Feasible means are
  untouched and the parity dumps are unchanged.

* `generate_design_analysis(file = )` wrote the script through a text-mode connection. On
  Windows the file arrived with CRLF line endings, which turns the SLURM part's first line
  into a `#!/bin/bash\r` shebang no cluster can execute, and on every platform it carried
  a doubled trailing newline. The script now reaches the disk byte for byte as returned,
  and the app's script and specification downloads write through the same binary path.

* The emitted SLURM wrapper could only run for its author. It hard-coded the author's
  cluster account and project paths, so any other user's submission failed at the
  scheduler while the surrounding instructions told them to save the parts under their own
  names. The wrapper now carries two placeholders marked `EDIT`, the `--account` directive
  and a writable `PROJECT_DIR`, and it invokes the `design_analysis.R` that sits beside it,
  wherever the user saved the pair.

* `precision_design()` documents all sixteen of its columns. The Monte Carlo standard
  errors and Wilson interval bounds were returned but missing from the reference page.

* A true effect of exactly zero returned `Inf`, and the package recommends that input.
  `design_conditions()` deliberately produces a null condition so a run can show how often
  it declares something when there is nothing to find. Putting it through `power_design()`
  returned `type_m = Inf`, straight into the app's display, while Type S silently
  degenerated to "the estimate is positive" because `true_effect > 0` is `FALSE` at zero.
  The Python twin raised `ZeroDivisionError` on the same input. Type S and Type M are now
  `NA` when the true effect is zero or unknown, in both engines, which is the guard
  `power_mixed()` already carried, applied to the other three sites.

* `response_variance()` reported a total that was not the sum of its parts. An undefined
  fixed component was laundered into zero while the result kept calling itself the sum. A
  single-row design now reports a fixed component of 0 and an honest total.

* The HPC precision-array job was broken and stale. It hand-listed the engine files and
  omitted `validate.R`, so it died on its first real line, and it used the indexed seed
  rule abandoned at 0.3, which correlates consecutive replicates and understates the job's
  own Monte Carlo error. It now sources the package wholesale, as the parity harness does,
  and uses the package's own replicate seeds, `qnorm(0.975)` and error reporting.

* A test now holds `print.pilotr_power()` to keeping its whole output on one stream. A
  header written with `message()` or through cli would land on the message stream while the
  table beneath it went to standard output, which knitr collects separately and renders as
  two boxes for one printed object. The method already used `cat()` throughout, and
  `test-print.R` keeps it there.

* Every vignette now turns console colour off and fixes the console width while it
  renders. pkgdown passes the calling terminal's colour support into its build
  subprocess, so a coloured message or error would otherwise reach the reader as escape
  sequences in the middle of the text.

## New in the generative core

* `varies_by = "observation"` draws a predictor once per row, for a quantity that varies trial by
  trial. `varies_by` is also validated now, because anything other than `"subject"` was previously
  read as item-level, so a predictor declared to vary by `"trial"` was silently given one value per
  item, and anyone who wrote that got wrong results.
* `dist = "uniform"`, with `min` and `max`. A uniform costs the same single draw as a normal, so
  it does not move the stream.
* `reliability` on a predictor simulates imperfect measurement. The latent value drives the
  linear predictor and any slope keyed on it, while the contaminated observed value goes into the
  data. No comparable package models unreliable predictors, and cross-level interactions are where
  unreliability bites hardest. Population moments are used rather than sample ones, since R's
  `mean()` and `sd()` accumulate in long double and 'Python's do not.
* The `exgaussian` family, the registered model family for reaction-time work, in brms's
  parameterisation so that a specification and the model fitted to it agree on what the intercept
  means. A shifted lognormal is not a substitute, because `model_data()` logs the response back and
  leaves a symmetric residual on the analysis scale.
* `response_variance()` decomposes the linear predictor's variance into the fixed part, each
  grouping factor's part, and the residual. Each grouping factor's component is exact for the
  realised design, averaging over the random-effect distribution analytically, because estimating
  it from the drawn effects of 30 subjects carries a sampling error of around a quarter of the
  component.

  A residual is reported for all eight families, not only the four carrying an explicit `sigma`,
  so the components are a complete decomposition everywhere and their ratios read as the design's
  intraclass correlations. For the link families it is the latent-scale distribution-specific
  variance (Nakagawa, Johnson and Schielzeth, 2017). Three of those four are derived from the
  process pilotr simulates and are exact for it: a `bernoulli` row is drawn as
  `1[u < invlogit(eta)]`, so its latent error is a standard logistic variate of variance `pi^2 / 3`,
  `ordinal` compares the same uniform against cumulative thresholds and inherits it, and for `beta`
  the identity `Var(logit(Y)) = trigamma(a) + trigamma(b)` is exact. `poisson` is the one
  approximation, since a count of zero has no logarithm. The trigamma form is used, and it is worth
  reading as an order of magnitude when counts are rare, where the published alternatives diverge
  sharply from it.
* `calibrate_response()` rescales a design to a target total variance, which is what lets a
  region of practical equivalence be stated in standard-deviation units and read the same way across
  designs. It now accounts for a residual it cannot move. A `bernoulli` or `ordinal` design carries
  a latent residual of about 3.29, so calibrating one to a total variance of 1 is impossible, and
  is refused with that number in the message. For `poisson` and `beta` the residual moves with the
  linear predictor, so the factor is solved numerically, and the whole search costs one extra
  simulation.

## New in the design-analysis layer

* `generate_design_analysis()` emits a Bayesian design analysis as a runnable script: a
  `brm()` call with `sample_prior = "yes"`, a Savage-Dickey Bayes factor, a highest-density
  interval against a region of practical equivalence, a three-way supported/null/inconclusive
  verdict, and a convergence gate that withholds every verdict when R-hat or the divergence rate
  fails it. Optionally it also emits a SLURM array wrapper and an aggregator. Emitting a script is
  what keeps this reachable from the browser build, where Stan cannot run at all.
* `power_mixed()` is no longer restricted to one within-unit factor and a crossed design, and no
  longer fits a formula written into the source. It takes `focal` and `formula`, derives the model
  from the specification through `model_formula()` and `model_data()`, and runs the same replicate
  loop as `precision_design()`. A design with two factors, with continuous predictors, or with a
  smaller random-effects structure than the maximal one was previously refused outright or analysed
  under a model it had not described.
* `sweep_spec()` runs an analysis over any addressed field of a specification.
  `power_curve_mixed()` and `precision_curve()` are now thin wrappers over it for sample size, which
  was previously the only axis reachable without writing the loop by hand. `design_conditions()`
  builds the coefficient sets for an effect-size sweep, including a shared all-zero condition.
* Every reported rate now carries its Monte Carlo standard error and a Wilson interval. The
  curve functions default to 60 replicates, at which a rate near 0.5 has a standard error of 0.065,
  which cannot support a claim anyone would want to make. The Wilson interval is used because the
  plain standard error collapses to zero at rates of 0 and 1, reading as certainty precisely where
  there is least of it.
* `precision_design()` and `power_mixed()` default `focal` to every coefficient in the
  specification.
* `print()` for a power result shows the fit accounting and each estimate beside its uncertainty.

* `solve_curve()` and `target_n()` solve a simulated design curve for the value that
  meets a target. The package computed the whole curve and then handed back the last and
  most consequential step, the number that goes into a preregistration, at replicate counts
  where neighbouring points are not significantly different. The plots drew a dashed line
  at 0.80 and left the reader to judge the crossing by eye. `target_n()` now takes the
  data frame `power_curve_mixed()` already returns and reports the sample size at which
  power reaches the target, with a confidence interval on it, rounded up to whole subjects.
  `solve_curve()` is the general form and reads `precision_curve()` and `sweep_spec()`
  output through the same column names, on any axis the sweep varied.

  The fit is a binomial probit regression of the decision rate against the swept value,
  weighted by the replicates behind each point and inverted by the delta method that
  `MASS::dose.p()` applies to a fitted `glm`. The probit follows from the design itself:
  under the normal approximation to a two-group comparison the probit of power is linear
  in the square root of the sample size. Checked against `stats::power.t.test()` over 36
  combinations of effect size, target power and grid shape, at 400 replicates a point, the
  solved sample size sat within 2.9% of the analytic answer on average and 6.9% at worst,
  against 3.4% and 8.9% for a logit fitted the same way.
  `tools/calibration/solve_curve_calibration.R` runs that comparison and writes every
  figure quoted for it to a file beside itself. Where the two-parameter model does not
  describe the curve, the interval is widened by the heterogeneity factor of probit
  analysis, reported as `dispersion`, since a curve the model fits badly should not report
  the narrow interval its replicate counts alone imply.

  Nothing extrapolates. A curve that does not reach the target within the sizes it swept is
  refused, with the range it did cover, and so is a fit that solves past the end of the
  sweep, a curve with no trend to invert, and a slope that cannot be told from zero. The
  app reports the solved size beside the curve it draws, and prints the refusal when the
  sweep does not settle the question.

* The Bayesian design-analysis record carries the rule that produced it, namely the
  decision thresholds, the convergence gate, a specification fingerprint and the package
  version, and the aggregator refuses to pool replicates that disagree on any of them.

## References

Barr, D. J., Levy, R., Scheepers, C., & Tily, H. J. (2013). Random effects structure for
confirmatory hypothesis testing: Keep it maximal. *Journal of Memory and Language, 68*(3), 255-278.
<doi:10.1016/j.jml.2012.11.001>

Nakagawa, S., Johnson, P. C. D. and Schielzeth, H. (2017). The coefficient of determination R2 and
intra-class correlation coefficient from generalized linear mixed-effects models revisited and
expanded. *Journal of the Royal Society Interface, 14*(134), 20170213.
<doi:10.1098/rsif.2017.0213>

# pilotr 0.2.1

Released 2026-07-23.

## Documentation

* A package-level help topic, `?pilotr`, introduces the workflow and groups
  every export by the stage it belongs to.
* The reference examples that cut the replicate count so they run quickly now
  say so, and each `n_sims` argument documents the Monte Carlo error and a
  usable recommendation for real planning.
* The `precision_curve()` example passes a region of practical equivalence
  clearly inside the effect it declares. It previously fell back to a default
  equal to that effect, the configuration the vignette warns against.

## Fixes

* `model_formula()` binds the global environment, so a formula prints without
  the environment tag that varied between builds and leaked into the rendered
  reference pages.

# pilotr 0.2.0

Released 2026-07-15.

## Design specifications

* `pilotr_example()` lists the design specifications shipped with the package, one per
  design family, and returns the path to each for `load_spec()`. They are the same JSON
  files that drive the Python package and the no-code app, and the Python package gains
  the same function.
* The eight per-family specifications now travel inside the package itself, so an
  installed copy can load them without a checkout of the repository.

## Documentation

* A new *Worked examples* article simulates every shipped specification.

# pilotr 0.1.0

Initial release.

## Design specifications and simulation

* `simulate_design()` generates a data set from a portable JSON design specification
  shared with the Python package of the same name, with `build_spec()` composing a
  specification from a flat list of design inputs, `load_spec()` reading one authored
  elsewhere (such as one downloaded from the no-code app) and `spec_json()` serialising
  it back to JSON.
* The response families are Gaussian, lognormal, shifted lognormal, Bernoulli, Poisson,
  ordinal and Beta.
* A shared cross-language random-number generator (`make_rng()`, with the AS 241 inverse
  normal in `as241()`) makes the simulated data bit-identical to the Python package's
  given the same specification and seed.
* A `per_subject` value below 1 or above the number of items is rejected with a clear
  error. The Python package raises on the same inputs.

## Power and precision analysis

* Simulation-based power and design analysis with `power_design()`, `power_mixed()` and
  `power_curve_mixed()`, reporting power alongside the Type S and Type M errors of
  Gelman and Carlin (2014). `power_curve_mixed()` returns an `n_converged` column
  alongside `power` and `type_m`, so the Monte Carlo standard error can be computed over
  the converged replicates.
* Precision and region-of-practical-equivalence (ROPE) design analysis with
  `precision_design()` and `precision_curve()`.
* Every analysis seeds the shared RNG as `seed + (replicate - 1)`, the same indexed-seed
  rule across all power and precision functions, the HPC array script and both languages.
* `power_mixed()` and `power_curve_mixed()` reject a specification without an item unit
  with a clear error instead of returning `NaN` power.
* All simulation-based power and precision analyses take a `workers` argument that
  spreads the Monte Carlo replicates across local cores with base R's `parallel`
  package. Because every replicate seeds the shared RNG from its own index, any worker
  count returns results identical to a serial run, and the sweep functions create their
  worker pool once and reuse it across all sample sizes.

## Modelling bridges and the no-code app

* `model_formula()` and `model_data()` derive the `lmer` formula and modelling data
  frame implied by a specification, and `brms_bridge()` derives a `brms` formula, family
  and priors, mapping every response family to its `brms` counterpart (including
  `Beta()`).
* `generate_r_script()` writes a self-contained, reproducible R script from a
  specification.
* A no-code Shiny application over the same specification, launched with `run_app()`.

## Metadata

* `citation("pilotr")` builds its version note from the package metadata, so it cannot
  drift from `DESCRIPTION`.
