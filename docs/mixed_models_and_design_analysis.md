# Advanced designs: continuous predictors, interactions, multilevel structure, and the Bayesian bridge

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
| Random slopes on continuous predictors | random `slopes` keyed by predictor, e.g. `(1 + x1 + x2 | subject)` |
| Additional grouping factors | any `random` entry beyond subject/item with `over` + `n` (e.g. subjects nested in clusters) |
| Positive/log outcomes | `lognormal` and `shifted_lognormal` families |
| Reaction times with a long right tail | the `exgaussian` family (`sigma`, `beta`), whose mean is the linear predictor, as in brms's `exgaussian()` |
| Bounded proportions | the `beta` family, with a mean and precision parameterisation |
| Partial crossing | a per-subject item count (`item.per_subject`), so each subject sees a sampled subset of items |
| Precision-based design analysis | `precision_design()` / `precision_curve()` report P(95% CI outside/inside a ROPE) and expected CI width, swept over *N* |
| Bayesian bridge | `brms_bridge()` returns a ready-to-fit `brms` formula and priors from the spec |

## Validation (`r/pilotr/examples/precision_design_analysis.R`, `r/pilotr/examples/run_demo.R`, `python/examples/run_demo.py`)

- Recovery. A maximal `lmer` fit to a large simulation recovers the full fixed-effect
  model, including main effects, interactions and the by-subject random slopes on continuous
  predictors (e.g. interaction 0.020→0.020, by-subject slope SDs 0.055/0.042 vs spec 0.06/0.04,
  residual 0.251 vs 0.25).
- Precision/ROPE design analysis. For a meaningful effect (β=0.10, outside the ROPE) the
  probability the 95% CI lands entirely outside the ROPE rises with *N*. For a negligible effect
  (β=0.02, inside the ROPE) the probability it lands entirely inside is high. This is the
  frequentist analogue of a Bayesian HDI-vs-ROPE decision.
- Additional grouping. A design with subjects nested in clusters generates a cluster random
  intercept (`nested_clusters.json`).
- Cross-language parity. The nested design and the continuous-predictor design
  (`reading_time_continuous.json`) give the same data in R and Python, to the last bit wherever
  the two share a maths library, as on Linux, where continuous integration checks every shipped
  example with the harness in `tools/parity/`. The cross-language article of the R package sets
  out what can differ on other platforms. Factor-only specs keep the original RNG stream.

## Recommended workflow

1. Build the design spec. You can point and click in the app (using the *advanced: paste a JSON
   spec* option for continuous-predictor designs), or edit a spec under `spec/examples/`. The
   same spec runs in R and Python.
2. Plan *N* by precision. `precision_curve(spec, focal, subject_ns, rope = 0.05)` sweeps the
   sample size and reports, per focal effect, P(meaningful) and P(equivalent) along with the
   expected CI width.
3. Solve the curve. `solve_curve(curve, target = 0.90)` fits the decision probability against
   the sample size and inverts the fit, returning the minimum analysable *N* for a determinate
   ROPE decision together with a confidence interval on it. `target_n(curve)` is the same thing
   for a power curve against a target power. Neither extrapolates: a curve that does not reach
   the target within the sizes swept is refused rather than extended.
4. Confirm in `brms`. `brms_bridge(spec)` returns the confirmatory model (family, fixed/random
   formula and weakly-informative priors, scaled to the response for the continuous families).
   Simulate one dataset with `model_data(spec, simulate_design(spec))`, which adds the contrast
   columns the formula reads, fit it in `brms` and check recovery, prior-predictive behaviour
   and the HDI-vs-ROPE decision before collecting data.

In this division of labour, pilotr serves as a design-and-planning engine, with `brms` (Stan)
providing the confirmatory Bayesian layer.

## Deferred (roadmap)

| Gap | Note |
|---|---|
| Full Bayesian fit in-toolkit (HDI, LKJ priors, posterior precision) | use the `brms_bridge()` output. The built-in precision analysis is a fast frequentist analogue |
| Correlated random effects and small-sample tests in Python | Python's `power_mixed` fits the specification's model with independent variance components and tests with a Wald z, which with few subjects or items rejects more readily than R's Satterthwaite tests. R's `power_mixed()` fits the same model with its correlations. Neither a joint test of several coefficients nor a design without random effects is planned in Python yet |
| Frequentist power for ordinal and Beta outcomes, and a GLM backend in Python | R's `power_mixed()` and `precision_design()` fit accuracy and counts with `glmer()`, but fit ordinal and Beta outcomes with a linear model on the response scale, warn and withhold Type M and the ROPE decisions. Python's `power_mixed` does the same for all four families. `generate_design_analysis()` analyses them on the link scale |
