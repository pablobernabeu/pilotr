# Advanced designs: continuous predictors, interactions, multilevel structure and the Bayesian bridge

Beyond simple factorial experiments, pilotr expresses observational and regression-style
mixed-effects designs. This note summarises the advanced capabilities, the features that are
validated, the features that are deferred and the recommended workflow for design and analysis
planning.

## Capabilities

| Capability | How |
|---|---|
| Continuous predictors | `predictors` block, with one value per subject, item or row (`varies_by` is `"subject"`, `"item"` or `"observation"`), drawn `N(mean, sd)` or, with `dist` set to `"uniform"`, `Uniform(min, max)` |
| Predictors measured with error | `reliability` on a predictor: the latent value drives the response, and the data carry an observed value that correlates `sqrt(reliability)` with it |
| Interactions | fixed-effect coefficient keys like `"x1:x2"` (product of columns) |
| Random slopes on continuous predictors | random `slopes` keyed by predictor, e.g. `(1 + x1 + x2 \| subject)` |
| Additional grouping factors | any `random` entry beyond subject/item with `over` + `n` (e.g. subjects nested in clusters) |
| Positive/log outcomes | `lognormal` and `shifted_lognormal` families |
| Reaction times with a long right tail | the `exgaussian` family (`sigma`, `beta`), whose mean is the linear predictor, as in brms's `exgaussian()` |
| Bounded proportions | the `beta` family, with a mean and precision parameterisation |
| Partial crossing | a per-subject item count (`units.item.per_subject`), so each subject sees a sampled subset of items |
| Precision-based design analysis | `precision_design()` / `precision_curve()` report P(95% CI outside/inside a ROPE) and expected CI width, swept over *N* |
| Bayesian bridge | `brms_bridge()` returns a ready-to-fit `brms` formula and priors from the spec |

## Validation (`r/pilotr/examples/precision_design_analysis.R`, `r/pilotr/examples/run_demo.R`, `python/examples/run_demo.py`)

- An `lmer` fit to a large simulation recovers the main effects and interactions of the full
  fixed-effect model and the SDs of the by-subject random slopes on continuous predictors (e.g.
  interaction 0.020→0.020, by-subject slope SDs 0.055/0.042 vs spec 0.06/0.04, residual 0.251
  vs 0.25). The fit uses the model that `model_formula()` derives from the spec. Its random
  effects are those of the generating model.
- In the precision/ROPE design analysis, the probability that the 95% CI lies entirely outside
  the ROPE rises with *N* for a meaningful effect (β=0.10, outside the ROPE). The probability that
  it lies entirely inside is high for a negligible effect (β=0.02, inside the ROPE). This is the
  frequentist analogue of a Bayesian HDI-vs-ROPE decision.
- A design with subjects nested in clusters (`nested_clusters.json`) generates a cluster random
  intercept.
- The nested design and the continuous-predictor design (`reading_time_continuous.json`) give the
  same data in R and Python, to the last bit wherever the two share a maths library, as on Linux.
  There, continuous integration checks every example spec with the harness in `tools/parity/`.
  The harness holds every case to the last bit except the unrounded dumps of `crossed_mixed_rt`,
  `reading_time_continuous` and `stress_exgaussian`, which `tools/parity/tolerance.json` allows
  to differ by up to 8 units in the last place, on Linux as on any other platform. The
  cross-language article of the R package sets out what can differ on other platforms.
  Factor-only specs keep the original RNG stream.

## Recommended workflow

1. Build the design spec. You can point and click in the app, using its advanced option to paste
   a JSON spec for continuous-predictor designs, or edit a spec under `spec/examples/`. The same
   spec runs in R and Python.
2. Plan *N* by precision. `precision_curve(spec, focal, subject_ns, rope = 0.05)` sweeps the
   sample size and reports, per focal effect, P(meaningful) and P(equivalent) along with the
   expected CI width.
3. Solve the curve. `solve_curve(curve, target = 0.90, effect = "SyntaxPC")` fits one focal
   effect's decision probability against the sample size and inverts the fit. It returns the
   minimum analysable *N* for a determinate ROPE decision together with a confidence interval on
   it. `effect` names the focal effect, here a coefficient of `reading_time_continuous.json`, and
   is required for a curve over several focal effects. `target_n(curve)` is the same thing for a
   power curve against a target power, and takes `effect` in the same way. Neither extrapolates,
   and both give an error for a curve that does not reach the target within the sizes swept.
4. Confirm in `brms`. `brms_bridge(spec)` returns the confirmatory model (family, fixed/random
   formula and weakly-informative priors, scaled to the response for the continuous families).
   Simulate one dataset with `model_data(spec, simulate_design(spec))`, which adds the contrast
   columns the formula reads. Fit it in `brms` (Stan) and check recovery, prior-predictive
   behaviour and the HDI-vs-ROPE decision before collecting data.

## Deferred (roadmap)

| Gap | Note |
|---|---|
| Full Bayesian fit in-toolkit (HDI, LKJ priors, posterior precision) | use the `brms_bridge()` output. The built-in precision analysis is a fast frequentist analogue |
| Correlated random effects and small-sample tests in Python | Python's `power_mixed` fits the spec's model with independent variance components and tests with a Wald z, which with few subjects or items rejects more readily than R's Satterthwaite tests. R's `power_mixed()` fits the same model with its correlations. Python's `power_mixed` has no joint test of several coefficients and raises an error for a design without random effects, both of which R's `power_mixed()` handles |
| Frequentist power for ordinal and Beta outcomes, and a GLM backend in Python | R's `power_mixed()` and `precision_design()` fit accuracy and counts with `glmer()`, but fit ordinal and Beta outcomes with a linear model on the response scale, warn and withhold Type M and the ROPE decisions. Python's `power_mixed` does the same for all four families. `generate_design_analysis()` analyses them on the link scale |
