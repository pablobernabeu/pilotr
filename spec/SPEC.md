# The `pilotr` design specification (v0.3)

A design specification is a single JSON object describing a data-generating process
(DGP) for an experiment. It is the contract shared by the web application, the R package, and the
Python package. Given the same spec and seed, every implementation must produce an
identical data set, within the boundary set out under **Scope of the guarantee** below.

The machine-readable form of this document is [`design.schema.json`](design.schema.json). Both
implementations enforce it through `validate_spec()`, together with the cross-field rules JSON
Schema cannot express.

## Versioning

`spec_version` is a `"major.minor"` string. A specification without it is a 0.2 specification,
which is what every specification written before the field existed is.

A specification that uses a feature introduced in 0.3 must declare `"0.3"` or later. This is not
bookkeeping: a 0.2 implementation reads such a specification differently and generates different
data while reporting success, which is the worst failure mode available to a reproducibility tool.
The 0.3 features are observation-level predictors, `dist`, `reliability`, the `exgaussian` family,
the `correlated` flag, and interaction random slopes. An implementation from 0.3 onwards refuses a
specification declaring a version newer than it understands, rather than reading part of it.

## Top-level fields

| Field | Type | Meaning |
|---|---|---|
| `spec_version` | string | Specification version, `"major.minor"`. Absent means 0.2. |
| `name` | string | Human label for the design. |
| `seed` | integer | Master seed (see RNG contract below). |
| `units` | object | Sampling units, e.g. `{"subject": {"n": 30}, "item": {"n": 24}}`. `item` is optional. Add `per_subject` to `item` (e.g. `{"n": 40, "per_subject": 12}`) for partial crossing, in which each subject sees a random subset of items. |
| `factors` | array | Experimental factors (categorical; see below). |
| `predictors` | array | Optional continuous predictors (see below). |
| `fixed` | object | Fixed effects: `intercept` + `coefficients` (map column → β). A coefficient key may be a single column or an `"a:b"` interaction (the product of columns a and b). |
| `random` | object | Random-effect structure by unit (`subject`, `item`). Empty `{}` ⇒ no random effects. |
| `response` | object | Outcome family + parameters (see below). |

A coefficient or slope key that names no existing column contributes zero, which silently removes
the term rather than failing. `validate_spec()` therefore refuses such a key: a design whose focal
effect is misspelled generates exactly the data of a null design and reports success.

### Factors

```json
{ "name": "condition",
  "levels": ["related", "unrelated"],
  "contrasts": { "cond": [-0.5, 0.5] },
  "vary_within": ["subject", "item"] }
```

* `contrasts` maps one or more contrast-column names to a numeric value per level
  (length = number of levels). Fixed coefficients and random slopes are keyed by these
  contrast-column names. This follows the convention used in `lme4` and in DeBruine and
  Barr (2021), where effects are coefficients on contrast-coded predictors.
* `vary_within`: the factor is crossed *within* the listed units (a within-unit factor),
  expanding each unit combination into one row per level.
* `between`: `"subject"` or `"item"`. The factor partitions that unit into equal blocks
  in level order (a between-unit factor that does not expand rows).

### Continuous predictors

```json
"predictors": [
  { "name": "SyntaxPC", "varies_by": "item", "mean": 0, "sd": 1 },
  { "name": "age", "varies_by": "subject", "mean": 0, "sd": 1 }
]
```

Each continuous predictor draws one value per unit and assigns it to all of that unit's rows. The
predictor name is a column usable in fixed `coefficients` (as a main effect or in an `"a:b"`
interaction) and in random-effect `slopes` (e.g. a by-subject random slope on an item-level
predictor, as in `(1 + SyntaxPC | subject)`). The defaults are `mean` 0 and `sd` 1.

`varies_by` is one of `"subject"`, `"item"` or `"observation"`. The last of these, new in 0.3,
draws one value per row, which is what a predictor varying trial by trial needs. Before 0.3 any
value other than `"subject"` was read as item-level, so a predictor declared to vary by `"trial"`
was silently given one value per item; it is now validated against the three names that exist.

`dist` (new in 0.3) selects the distribution, either `"normal"`, which uses `mean` and `sd`, or
`"uniform"`, which uses `min` and `max`. A uniform draw consumes exactly as much of the random
stream as a normal one, since a normal is produced by transforming a single uniform, so switching
between them does not move the stream.

#### Reliability

`reliability` (new in 0.3) simulates a predictor measured with error:

```json
{ "name": "z_reading", "varies_by": "subject", "sd": 1, "reliability": 0.8 }
```

The **latent** value drives the linear predictor and any random slope keyed on the predictor, while
the **observed**, contaminated value is what appears in the returned data, which is what an analyst
would actually have measured. Writing `ρ` for the reliability and using the predictor's population
mean and standard deviation,

```
observed = mean + (true − mean + sd·sqrt((1 − ρ)/ρ)·z) · sqrt(ρ)
```

so the observed variable has the same variance as the latent one and correlates `sqrt(ρ)` with it.
Reliability in the classical sense is that squared correlation, which is why the field is `ρ`.

The attenuation is `sqrt(ρ)` rather than the `ρ` of the textbook regression-dilution result because
both variables are placed on the same variance here; standardising the observed variable back to
the latent one's variance absorbs the `1/sqrt(ρ)` factor that result carries.

The moments used are the **population** ones, not the sample mean and standard deviation of the
values drawn. R's `mean()` and `sd()` accumulate in long double and Python's do not, so
standardising against the sample would reintroduce a cross-language divergence.

One further normal is drawn per value, and only when `reliability` is present and below 1, so a
specification that does not use it keeps the original stream. No comparable package models
unreliable predictors, and cross-level interactions are where unreliability bites hardest.

### Random effects (per unit)

```json
"subject": {
  "intercept_sd": 0.12,
  "slopes": { "cond": 0.04 },
  "correlations": { "intercept,cond": 0.2 }
}
```

The random-effect column order is `["intercept", <slopes in listed order>]`. A
covariance matrix `Σ = D · R · D` is formed from the SDs `D` and the correlation matrix `R`,
which has a unit diagonal and off-diagonals taken from `correlations`, keyed `"a,b"` (a tilde
separator, `"a~b"`, is also accepted). Each element is computed as `(sd_i · sd_j) · r_ij`; the
bracketing is part of the contract, because floating-point multiplication is not associative and
the alternative grouping lands on a different double for some inputs. Per unit, a vector
`b = L z` is drawn, where `L` is the lower Cholesky factor of `Σ` and `z` are iid standard
normals. The unit's contribution to a row's linear predictor is
`b[intercept] + Σ_k b[slope_k] · (design value of slope_k for that row)`.

Slope keys follow exactly the same rule as fixed coefficients: a contrast column, a continuous
predictor, or an `"a:b"` interaction between them. Before 0.3 an interaction slope was accepted,
sized into the covariance, drawn from the stream, and then silently discarded, so the emitted
analysis model contained a term the generative process did not, which inflates power in the
direction Barr et al. (2013) warn about.

`Σ` must be positive definite. One that is not is an error naming the grouping factor and the
random-effect column at which the factorisation failed. Clamping the failing pivot at zero and
continuing, as earlier versions did, produced random effects whose standard deviations were several
times the requested ones without reporting anything. A standard deviation of exactly zero is a
different matter and remains valid, since it is how a term is held fixed while the rest of the
structure is kept intact.

`correlated` (new in 0.3) says whether a group's random effects are correlated, defaulting to
whether `correlations` is supplied. It decides whether an emitted `lmer` or `brms` formula uses a
single or a double bar. Earlier versions always emitted a single bar and an LKJ prior, so a design
whose slopes were uncorrelated by construction was nonetheless analysed as though a correlation were
there to estimate. Setting `correlated` to `false` while also supplying `correlations` is
contradictory and is refused.

### Additional grouping factors

Any `random` entry whose name is not `subject` or `item` is an extra grouping factor. It
adds `over` (the unit it groups, either `"subject"` or `"item"`) and `n` (the number of groups).
The units are assigned to groups in equal blocks. For example, subjects nested in clusters:

```json
"site": { "over": "subject", "n": 12, "intercept_sd": 0.5, "slopes": { ... } }
```

Each group draws a random-effect vector (intercept + any slopes) applied to all rows of the
units in that group, and the simulated data gains a column with the group id. Useful for
hierarchical designs (e.g. participants within sites, schools, or languages).

### Response families

| `family` | Parameters | Generation |
|---|---|---|
| `gaussian` | `sigma` | `y = η + σ·z` |
| `shifted_lognormal` | `sigma`, `shift` | `y = shift + exp(η + σ·z)` (reaction times) |
| `lognormal` | `sigma` | `y = exp(η + σ·z)` (positive outcomes, e.g. reading time per word) |
| `exgaussian` | `sigma`, `beta` | `y = η + σ·z − β·(log(u) + 1)` (reaction times; new in 0.3) |
| `bernoulli` | — | `p = invlogit(η)`, `y = 1[u < p]` (accuracy; logit link) |
| `poisson` | — | `λ = exp(η)`, `y =` inverse-CDF Poisson (counts; log link) |
| `ordinal` | `thresholds` (K−1 cut-points) | cumulative-logit: `P(Y≤k) = invlogit(θ_k − η)` (Likert) |
| `beta` | `phi` (precision) | `μ = invlogit(η)`, `y ~ Beta(μ·φ, (1−μ)·φ)` (proportions in (0,1)) |

The ex-Gaussian is a normal plus an exponential, mean-centred by subtracting the exponential's own
mean, so that `η` remains the mean of the response. That is brms's `exgaussian(mu, sigma, beta)`
parameterisation, in which `mu` is the mean, so a specification and the model fitted to it agree on
what the intercept means. `−log(u)` is a unit exponential, hence the two draws per row. A shifted
lognormal is not a substitute, because `model_data()` logs the response back and leaves a symmetric
residual on the analysis scale.

`η` (the linear predictor for a row) = `intercept + Σ β_key · value(key)`, where a key is a
contrast column, a continuous predictor, or an `"a:b"` interaction (the product of the named
columns), `+` subject random part `+` item random part `+` the random parts of any additional
grouping factors. `name` sets the output column name. An optional `round` sets the decimal
rounding of the response, and applies only to the families whose outcome is continuous, the others
being integers already.

## Accumulation order (identical across all implementations)

Floating-point addition is not associative, so `(a + b) + c` and `a + (b + c)` can land on different
doubles. Summing the terms of `η` in a different order therefore produces a different result, and a
guarantee of identical data has to fix the order as firmly as it fixes the draw order. Measured over
200,000 rows at realistic coefficient magnitudes, two orderings that differ only in their bracketing
disagreed on 63.8% of rows.

`η` accumulates as a strict left fold, one term at a time, in this order:

1. `η ← intercept`
2. for each fixed coefficient, in the order the `coefficients` object lists them:
   `η ← η + β_key · value(key)`
3. if the subject group exists: `η ← η + b[intercept]`, then for each subject slope in listed
   order, `η ← η + b[slope_k] · value(slope_k)`
4. the same for the item group, if it exists
5. the same for each additional grouping factor, in the order the `random` entries are listed

An interaction value is itself a left fold: `v ← 1`, then `v ← v · value(part)` for each part of the
key in written order.

Neither language's built-in summation may be used for any of this. Base R's `sum()` accumulates in
80-bit long double on x86, and CPython's `sum()` has applied Neumaier compensation to floats since
version 3.12. Both are more accurate than a plain double fold, but they are more accurate in
different ways, so an inner product of three terms or more can land on different doubles in the two
ports. Every inner product, including those inside the Cholesky factorisation and the matrix-vector
product, is written as an explicit double fold.

For the same reason, integer powers are written as repeated multiplication rather than with `^` or
`**`. R special-cases small integer exponents while Python calls the library `pow()`; measured over
200,000 draws in the Gamma sampler's range, the two disagreed on a third of inputs by up to 6 ulp,
and because that value decides a rejection step, the disagreement changed how many draws were
consumed.

## RNG contract (identical across all implementations)

**Uniform generator.** L'Ecuyer (1988) combined LCG:

```
s1 ← (40014 · s1) mod 2147483563
s2 ← (40692 · s2) mod 2147483399
d  ← s1 − s2 ;  if d < 1 then d ← d + 2147483562
u  ← d / 2147483563            # u ∈ (0, 1)
```

All products stay below 2^53, so the arithmetic is exact in IEEE-754 doubles and in
Python integers alike. The seeding rule is `s1 ← 1 + (|seed| mod 2147483562)` and
`s2 ← 1 + ((40692 · s1) mod 2147483398)`, after which 10 warm-up draws are discarded.

**Normal deviates.** Wichura (1988) Algorithm AS 241 applied to `u`, the algorithm R's
`qnorm` uses. Deviates therefore agree to full double precision.

**Draw order (must be identical everywhere):**

0. If `units.item.per_subject` is set (partial crossing): for each subject `s = 1..S`, in
   row-build order, sample that subject's item subset by a partial Fisher–Yates shuffle,
   consuming one uniform per sampled item (`per_subject` uniforms per subject). These are
   the first RNG draws. (Skipped entirely under full crossing, so fully crossed specs keep
   the original stream.)
1. For each continuous predictor (in listed order): for each of its units `u = 1..N`, draw one
   deviate, `N(mean, sd)` by default or `Uniform(min, max)` when `dist` is `"uniform"`, and then,
   only when `reliability` is present and below 1, one further standard normal for that value's
   measurement error. `N` is the number of subjects, of items, or of rows, according to
   `varies_by`. (Skipped entirely when there is no `predictors` block, so factor-only specs keep
   the original stream. A uniform costs the same one draw as a normal, and a `reliability` of 1 or
   absent costs nothing, so neither moves the stream either.)
2. For each subject `s = 1..S`: draw `q_subject` standard normals (intercept, then each
   slope in listed order); set `b_subject[s] = L_subject · z`.
3. For each item `t = 1..I` (if items exist): draw `q_item` standard normals; set
   `b_item[t] = L_item · z`.
4. For each additional grouping factor (in the order the `random` entries are listed): for
   each group `g = 0..K−1`, draw `q_group` standard normals; set `b_group[g] = L_group · z`.
5. Iterate observations in **canonical row order** and draw the response. Most families
   consume exactly one deviate per row (normal for gaussian/lognormal/shifted_lognormal;
   uniform for bernoulli/poisson/ordinal). The ex-Gaussian consumes exactly two, a normal then a
   uniform, in that order. The beta family instead consumes a variable,
   data-dependent number of draws per row: two Gamma variates through the Marsaglia–Tsang
   rejection sampler, each consuming normal–uniform pairs until acceptance (with one extra
   uniform per Gamma variate whose shape is below 1).

**Extending the draw order.** A new feature must consume **zero** draws when it is not used, so
that every specification written before it stays bit-identical. Each of the 0.3 additions is built
that way: an observation-level predictor only appears when declared, a uniform draw costs the same
as the normal it replaces, a `reliability` of 1 or absent draws nothing, and a family branch is
reached only when that family is selected.

**Replicate seeds.** The power and precision loops derive their per-replicate seeds from the
specification's seed by drawing from the shared generator, skipping any duplicate, rather than by
adding the replicate index. Consecutive seeds are not independent streams here: seeding sets `s1` to
`1 + (seed mod 2147483562)` and `s2` from `s1`, with only ten warm-up draws discarded, so the first
draw of replicate `i` correlated 0.95 with that of replicate `i + 1` under the old rule. An
arithmetic scramble does not help, since the seeding rule is linear in the seed. This changed every
number pilotr produced before 0.3.

## Scope of the guarantee

Identical data means bit-identical, and holds exactly for:

* the `gaussian` family, and any design whose response path applies no transcendental function to
  the linear predictor;
* **any** family when `response.round` is set, since rounding quantises away a last-bit difference.

For `lognormal`, `shifted_lognormal`, `exgaussian`, `bernoulli`, `poisson`, `ordinal` and `beta`
without `round`, results may differ in the last unit in the last place. IEEE-754 requires correct
rounding for addition, subtraction, multiplication, division and square root, but **not** for
`exp()` or `log()`, and the R and Python builds on a given platform need not share a maths library.
Measured over 200,000 arguments in the log-reaction-time range, R and CPython `exp()` disagreed on
0.44% of them by up to 6 ulp, and `log()` on 0.12% by up to 1 ulp.

This is demonstrable rather than assumed. Taking a shifted-lognormal design and switching only its
family to `gaussian`, so that the seed, the random-effect structure, the linear predictor and the
entire draw sequence are unchanged, gives bit-identical output in both languages, while the
lognormal original differs in a handful of rows. `exp()` is the only remaining difference.

For the discrete families the practical consequence is different in kind. A last-bit difference in
`exp()` usually changes nothing, because the outcome is an integer decided by a comparison; but when
the comparison sits exactly on a threshold, the outcome moves by a whole category. This is rare and
it is not impossible, so a design analysis that has to be reproducible to the last observation should
either set `round` or stay with `gaussian`.

The stricter guarantee within one language is unconditional: the same implementation, specification
and seed always produce the same data.

**Canonical row order:** nested loops, outermost first,
`for s in 1..S: for t in 1..I: for (each within-factor level-combination, factors in
listed order, levels in listed order): emit row`. Between-unit factors assign a level to
each unit by equal blocks in level order and do not expand rows.

## References

* L'Ecuyer, P. (1988). Efficient and portable combined random number generators.
  *Communications of the ACM, 31*(6), 742–751. https://doi.org/10.1145/62959.62969
* Wichura, M. J. (1988). Algorithm AS 241: The percentage points of the normal
  distribution. *Applied Statistics, 37*(3), 477–484. https://doi.org/10.2307/2347330
* DeBruine, L. M., & Barr, D. J. (2021). Understanding mixed-effects models through data
  simulation. *Advances in Methods and Practices in Psychological Science, 4*(1).
  https://doi.org/10.1177/2515245920965119
* Barr, D. J., Levy, R., Scheepers, C., & Tily, H. J. (2013). Random effects structure for
  confirmatory hypothesis testing: Keep it maximal. *Journal of Memory and Language, 68*(3),
  255–278. https://doi.org/10.1016/j.jml.2012.11.001
* Matuschek, H., Kliegl, R., Vasishth, S., Baayen, H., & Bates, D. (2017). Balancing Type I error
  and power in linear mixed models. *Journal of Memory and Language, 94*, 305–315.
  https://doi.org/10.1016/j.jml.2017.01.001
* Neumaier, A. (1974). Rundungsfehleranalyse einiger Verfahren zur Summation endlicher Summen.
  *Zeitschrift für Angewandte Mathematik und Mechanik, 54*(1), 39–51.
  https://doi.org/10.1002/zamm.19740540106
