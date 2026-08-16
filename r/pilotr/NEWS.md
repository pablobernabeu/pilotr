# pilotr (development version)

* **A Poisson mean beyond the sampler's reach was returned as the iteration cap.** The
  inverse-CDF walk starts from `exp(-mean)`, which underflows to exactly zero once the
  mean passes about 746 — a poisson intercept of 7 already implies a mean of exp(7),
  about 1097 — after which every simulated count came back as the cap of 1000000 while
  reporting success. Both engines now refuse such a mean, naming `exp(eta)` and the
  offending value, with message text character-for-character identical across the twins.
  Feasible means are untouched and the parity dumps are unchanged.
* **`generate_design_analysis(file = )` wrote the script through a text-mode
  connection.** On Windows the file arrived with CRLF line endings, which turns the
  SLURM part's first line into a `#!/bin/bash\r` shebang no cluster can execute, and on
  every platform it carried a doubled trailing newline. The script now reaches the disk
  byte for byte as returned, and the app's script and specification downloads write
  through the same binary path.
* **The emitted SLURM wrapper could only run for its author.** It hard-coded the
  author's cluster account and project paths, so any other user's submission failed at
  the scheduler while the surrounding instructions told them to save the parts under
  their own names. The wrapper now carries two placeholders marked `EDIT` — the
  `--account` directive and a writable `PROJECT_DIR` — and it invokes the
  `design_analysis.R` saved next to it rather than a path on the author's cluster.
* `precision_design()` documents all sixteen of its columns. The Monte Carlo standard
  errors and Wilson interval bounds were returned but missing from the reference page.
* `spec_from_model()` gains test coverage: the recovered specification's units,
  `between`/`vary_within` placement, interaction keys read back off product columns,
  and random-effect estimates are checked against the design that generated the pilot
  data, alongside the refusal paths for models the function cannot read.
* **A true effect of exactly zero returned `Inf`, and the package recommends that input.**
  `design_conditions()` deliberately produces a null condition so a run can show how often
  it declares something when there is nothing to find. Putting it through `power_design()`
  returned `type_m = Inf`, straight into the app's display, while Type S silently
  degenerated to "the estimate is positive" because `true_effect > 0` is `FALSE` at zero.
  The Python twin raised `ZeroDivisionError` on the same input. Type S and Type M are now
  `NA` when the true effect is zero or unknown, in both engines — the guard `power_mixed()`
  already carried, applied to the other three sites.
* **`response_variance()` reported a total that was not the sum of its parts.** An
  undefined fixed component was laundered into zero while the result kept calling itself
  the sum. A single-row design now reports a fixed component of 0 and an honest total.
* **The declared minimum R version is now stated.** `Depends: R (>= 4.0.0)` — 4.0.0 is the
  release whose `round()` the cross-language claim assumes, so the floor is load-bearing
  rather than conventional. A new CI job checks the declared Suggests floors, which the
  matrix (release, devel, oldrel-1) sat well above and so never exercised.
* **The Bayesian design-analysis record carries the rule that produced it** — the decision
  thresholds, the convergence gate, a specification fingerprint and the package version —
  and the aggregator refuses to pool replicates that disagree on any of them.
* **The HPC precision-array job was broken and stale.** It hand-listed the engine files and
  omitted `validate.R`, so it died on its first real line; it also used the indexed seed
  rule abandoned at 0.3, which correlates consecutive replicates and understates the job's
  own Monte Carlo error. It now sources the package wholesale, as the parity harness does,
  and uses the package's own replicate seeds, `qnorm(0.975)` and error reporting.
* Four small divergences between the two engines are closed: a non-whole seed truncates
  identically, a whole `spec_version` reads the same however written, a non-object unit is
  reported rather than crashing with a base error, and the message text of each is now
  character-for-character identical across the twins.

* **A test that the packaged example specifications match the repository's own.**
  `spec/examples/*.json` is canonical and both packages carry a mirror so an installed
  copy can reach it, but nothing enforced the mirror. A load-and-simulate test cannot: a
  stale packaged copy still loads and simulates perfectly well, it simply describes a
  different design from the one the repository documents. The new `test-examples.R`
  compares the bytes, and is twinned with the Python suite's `test_examples.py`. It skips
  when the package is checked in isolation from the repository.

* **A test that `print.pilotr_power()` keeps its whole output on one stream.** A header
  written with `message()` or through cli would land on the message stream while the table
  beneath it went to standard output, which knitr collects separately and renders as two
  boxes for one printed object. The method already used `cat()` throughout;
  `test-print.R` now holds it there.
* Every vignette now turns console colour off and fixes the console width while it
  renders. pkgdown passes the calling terminal's colour support into its build
  subprocess, so a coloured message or error would otherwise reach the reader as escape
  sequences in the middle of the text.

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
last place, on 0.09% and 0.05% of cells. That residue is the libm limit described immediately below
rather than a defect, and `tools/parity/tolerance.json` records which cases are allowed it and why.

* The scope of the guarantee is now stated honestly in `spec/SPEC.md`. It is exact for
  `gaussian`, for any design applying no transcendental function to the linear predictor, and for
  any family with `response.round` set. The families that apply `exp()` or `log()` may differ in the
  last unit in the last place, because IEEE-754 does not require correct rounding for those
  functions and the two builds need not share a maths library. Measured rates are given, and the
  attribution is demonstrated rather than asserted. The same design switched to `gaussian`, with an
  identical seed, structure and draw sequence, is bit-identical.

## Validation and versioning

* `validate_spec()` checks a specification against the schema and against the cross-field rules
  the schema cannot express, and `load_spec()` now calls it by default. A strict draft-07 schema had
  shipped since 0.1 with no code path consulting it.
* Validation exists because several ways of getting a specification wrong produced plausible data
  rather than an error. A mistyped coefficient key resolved to no column and so silently set that
  effect to zero, which generates exactly the data of a null design and reports success. A
  response parameter left over from another family was ignored. Both are now refused.
* `spec_version` negotiates the format. A specification with no such field is a 0.2
  specification. One using a 0.3 feature must declare 0.3, because a 0.2 implementation reads it
  differently and generates different data without complaint. One declaring a version newer than
  the implementation understands is refused rather than partly read.
* `simulate_design()` gains `validate`, defaulting to `TRUE`. The replicate loops validate once and
  then skip it, so a sweep pays the cost once rather than once per replicate.

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
  `vary_within` or a single ordinal threshold stays an array instead of collapsing to a scalar.
* `generate_r_script()` embedded the specification through `deparse()`, which prints 15
  significant digits and so does not round-trip. It now emits numbers at full precision, which
  matters because the point of the script is bit-for-bit reproduction.
* `model_formula()` and `brms_bridge()` emitted correlated random effects unconditionally, along
  with an LKJ prior, while the generative process only correlates them when `correlations` is
  supplied. They now follow the new `correlated` flag and emit a double bar otherwise. A group with
  no slopes keeps a single bar, since `lme4` cannot parse `(1 || g)`.
* `model_data()` did not create the product column an interaction random slope needs, so the
  emitted formula referred to a variable the modelling data lacked. It now covers the union of the
  fixed-coefficient and random-slope keys.
* A focal effect that never appears in any fit now warns instead of returning decision
  proportions of zero, which read as 'this design can decide nothing' when the cause was a name that
  did not match the model.
* A replicate loop in which no fit succeeds now reports why, passing on the fitter's own
  message. An unidentifiable random-effects structure previously produced a silent result of `NA`.
* A correlation naming a random-effect term that does not exist is an error rather than a subscript
  failure.

## New in the generative core

* `varies_by = "observation"` draws a predictor once per row, for a quantity that varies trial by
  trial. `varies_by` is also validated now, because anything other than `"subject"` was previously
  read as item-level, so a predictor declared to vary by `"trial"` was silently given one value per
  item.
  Anyone who wrote that had wrong results.
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
  realised design, averaging over the random-effect distribution analytically rather than drawing
  from it, because estimating it from the drawn effects of 30 subjects carries a sampling error of
  around a quarter of the component.

  A residual is reported for all eight families, not only the four carrying an explicit `sigma`,
  so the components are a complete decomposition everywhere and their ratios read as the design's
  intraclass correlations. For the link families it is the latent-scale distribution-specific
  variance (Nakagawa, Johnson and Schielzeth, 2017). Three of those four are exact for the process
  pilotr actually simulates rather than borrowed conventions: a `bernoulli` row is drawn as
  `1[u < invlogit(eta)]`, so its latent error is a standard logistic variate of variance `pi^2 / 3`,
  `ordinal` compares the same uniform against cumulative thresholds and inherits it, and for `beta`
  the identity `Var(logit(Y)) = trigamma(a) + trigamma(b)` is exact. `poisson` is the one
  approximation, since a count of zero has no logarithm. The trigamma form is used, and it is worth
  reading as an order of magnitude when counts are rare, where the published alternatives diverge
  sharply from it.
* `calibrate_response()` rescales a design to a target total variance, which is what lets a
  region of practical equivalence be stated in standard-deviation units and read the same way across
  designs. It now accounts for a residual it cannot move. A `bernoulli` or `ordinal` design carries
  a latent residual of about 3.29, so calibrating one to a total variance of 1 is impossible and is
  refused with that number rather than silently missing the target. For `poisson` and `beta` the
  residual moves with the linear predictor, so the factor is solved numerically at the cost of one
  extra simulation rather than one per candidate.

## New in the design-analysis layer

* `generate_design_analysis()` emits a runnable Bayesian design analysis rather than running
  one: a `brm()` call with `sample_prior = "yes"`, a Savage-Dickey Bayes factor, a highest-density
  interval against a region of practical equivalence, a three-way supported/null/inconclusive
  verdict, and a convergence gate that withholds every verdict when R-hat or the divergence rate
  fails it. Optionally it also emits a SLURM array wrapper and an aggregator. Emitting a script
  rather than fitting a model is what keeps this reachable from the browser build, where Stan cannot
  run at all.
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
