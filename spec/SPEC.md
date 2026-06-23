# The `pilotr` design specification (v0.1)

A design specification is a single JSON object describing a data-generating process
(DGP) for an experiment. It is the contract shared by the web application, the R package, and the
Python package. Given the same spec and seed, every implementation must produce an
identical data set.

## Top-level fields

| Field | Type | Meaning |
|---|---|---|
| `name` | string | Human label for the design. |
| `seed` | integer | Master seed (see RNG contract below). |
| `units` | object | Sampling units, e.g. `{"subject": {"n": 30}, "item": {"n": 24}}`. `item` is optional. Add `per_subject` to `item` (e.g. `{"n": 40, "per_subject": 12}`) for partial crossing, in which each subject sees a random subset of items. |
| `factors` | array | Experimental factors (categorical; see below). |
| `predictors` | array | Optional continuous predictors (see below). |
| `fixed` | object | Fixed effects: `intercept` + `coefficients` (map column → β). A coefficient key may be a single column or an `"a:b"` interaction (the product of columns a and b). |
| `random` | object | Random-effect structure by unit (`subject`, `item`). Empty `{}` ⇒ no random effects. |
| `response` | object | Outcome family + parameters (see below). |

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

Each continuous predictor draws one value per unit (`varies_by` = `"subject"` or
`"item"`) from `N(mean, sd)` and assigns it to all of that unit's rows. The predictor name is
a column usable in fixed `coefficients` (as a main effect or in an `"a:b"` interaction) and in
random-effect `slopes` (e.g. a by-subject random slope on an item-level predictor, as in
`(1 + SyntaxPC | subject)`). The defaults are `mean` 0 and `sd` 1.

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
which has a unit diagonal and off-diagonals taken from `correlations`, keyed `"a,b"`. Per unit, a vector
`b = L z` is drawn, where `L` is the lower Cholesky factor of `Σ` and `z` are iid standard
normals. The unit's contribution to a row's linear predictor is
`b[intercept] + Σ_k b[slope_k] · (contrast value of slope_k for that row)`.

Slopes may be keyed by a contrast column or by a continuous predictor (e.g. a by-subject
random slope on an item-level predictor, `(1 + SyntaxPC | subject)`).

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
| `bernoulli` | — | `p = invlogit(η)`, `y = 1[u < p]` (accuracy; logit link) |
| `poisson` | — | `λ = exp(η)`, `y =` inverse-CDF Poisson (counts; log link) |
| `ordinal` | `thresholds` (K−1 cut-points) | cumulative-logit: `P(Y≤k) = invlogit(θ_k − η)` (Likert) |
| `beta` | `phi` (precision) | `μ = invlogit(η)`, `y ~ Beta(μ·φ, (1−μ)·φ)` (proportions in (0,1)) |

`η` (the linear predictor for a row) = `intercept + Σ β_col · contrast_col +`
subject random part `+` item random part. `name` sets the output column name. An optional
`round` sets the decimal rounding of the response.

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

1. For each continuous predictor (in listed order): for each of its units `u = 1..N`, draw
   one `N(mean, sd)` deviate. (Skipped entirely when there is no `predictors` block, so
   factor-only specs keep the original stream.)
2. For each subject `s = 1..S`: draw `q_subject` standard normals (intercept, then each
   slope in listed order); set `b_subject[s] = L_subject · z`.
3. For each item `t = 1..I` (if items exist): draw `q_item` standard normals; set
   `b_item[t] = L_item · z`.
4. Iterate observations in **canonical row order** and draw exactly one response deviate
   (normal for gaussian/lognormal/shifted_lognormal; uniform for bernoulli/poisson/ordinal)
   per row.

**Canonical row order:** nested loops, outermost first,
`for s in 1..S: for t in 1..I: for (each within-factor level-combination, factors in
listed order, levels in listed order): emit row`. Between-unit factors assign a level to
each unit by equal blocks in level order and do not expand rows.

## References

* L'Ecuyer, P. (1988). Efficient and portable combined random number generators.
  *Communications of the ACM, 31*(6), 742–751.
* Wichura, M. J. (1988). Algorithm AS 241: The percentage points of the normal
  distribution. *Applied Statistics, 37*(3), 477–484.
* DeBruine, L. M., & Barr, D. J. (2021). Understanding mixed-effects models through data
  simulation. *Advances in Methods and Practices in Psychological Science, 4*(1).
