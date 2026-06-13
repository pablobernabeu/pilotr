# Advanced designs — continuous predictors, interactions, multilevel structure, and the Bayesian bridge

Beyond simple factorial experiments, simdgp expresses observational and regression-style
mixed-effects designs. This note summarises the advanced capabilities, what is validated, what
is deferred, and the recommended workflow for design and analysis planning.

## Capabilities

| Capability | How |
|---|---|
| **Continuous predictors** | `predictors` block — subject- or item-level, drawn `N(mean, sd)` |
| **Interactions** | fixed-effect coefficient keys like `"x1:x2"` (product of columns) |
| **Random slopes on continuous predictors** | random `slopes` keyed by predictor, e.g. `(1 + x1 + x2 | subject)` |
| **Additional grouping factors** | any `random` entry beyond subject/item with `over` + `n` (e.g. subjects nested in clusters) |
| **Positive/log outcomes** | `lognormal` and `shifted_lognormal` families |
| **Precision-based design analysis** | `precision_design()` / `precision_curve()` — P(95% CI outside/inside a ROPE) + expected CI width, swept over N |
| **Bayesian bridge** | `brms_bridge()` emits a ready-to-fit `brms` formula + priors from the spec |

## Validation (`examples/precision_design_analysis.R`, `examples/run_demo.*`)

- **Recovery**: a maximal `lmer` fit to a large simulation recovers the full fixed-effect
  model — main effects, interactions, and the by-subject random slopes on continuous
  predictors (e.g. interaction 0.020→0.020; by-subject slope SDs 0.055/0.042 vs spec 0.06/0.04;
  residual 0.251 vs 0.25).
- **Precision/ROPE design analysis**: for a meaningful effect (β=0.10, outside the ROPE) the
  probability the 95% CI lands entirely outside the ROPE rises with N; for a negligible effect
  (β=0.02, inside the ROPE) the probability it lands entirely inside is high — the frequentist
  analog of a Bayesian HDI-vs-ROPE decision.
- **Additional grouping**: a design with subjects nested in clusters generates a cluster random
  intercept and is bit-identical across R and Python (`nested_clusters.json`; parity diff 0.0).
- **Cross-language parity**: the continuous-predictor design is bit-identical R≡Python
  (~5e-15 over 4,000 rows × 7 columns); factor-only specs keep the original RNG stream.

## Recommended workflow

1. **Build the design spec** — point-and-click in the app (use *advanced: paste a JSON spec*
   for continuous-predictor designs), or edit a spec under `spec/examples/`. The same spec runs
   in R and Python.
2. **Plan N by precision** — `precision_curve(spec, formula, focal, subject_ns, rope)` sweeps the
   sample size and reports, per focal effect, P(meaningful) / P(equivalent) and the expected CI
   width, so you can read off the minimum analysable N for a determinate ROPE decision.
3. **Confirm in `brms`** — `brms_bridge(spec)` emits the confirmatory model (family + fixed/random
   formula + weakly-informative priors). Simulate one dataset, fit it in `brms`, and check
   recovery, prior-predictive behaviour, and the HDI-vs-ROPE decision before collecting data.

This makes simdgp a **design-and-planning engine** with `brms` (Stan) as the **confirmatory
Bayesian layer** — an honest division of labour rather than re-implementing Stan.

## Deferred (roadmap)

| Gap | Note |
|---|---|
| Full Bayesian fit in-toolkit (HDI, LKJ priors, posterior precision) | use the `brms_bridge()` output; the built-in precision analysis is a fast frequentist analog |
| Beta family (bounded proportions) | not yet; Bernoulli/Poisson/ordinal are supported |
| Unbalanced / partial crossing (each subject sees a subset of items) | not yet; current designs fully cross subject × item |
| Crossed-random-slope power in Python | `statsmodels` overstates slope variance (conservative); R/`lme4` is the reference |
