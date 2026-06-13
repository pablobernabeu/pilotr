# The `simdgp` design specification (v0.1)

A design specification is a single JSON object describing a **data-generating process**
(DGP) for an experiment. It is the contract shared by the web app, the R package, and the
Python package. Given the same spec and seed, every implementation must produce the
**identical** data set.

## Top-level fields

| Field | Type | Meaning |
|---|---|---|
| `name` | string | Human label for the design. |
| `seed` | integer | Master seed (see RNG contract below). |
| `units` | object | Sampling units, e.g. `{"subject": {"n": 30}, "item": {"n": 24}}`. `item` is optional. |
| `factors` | array | Experimental factors (see below). |
| `fixed` | object | Fixed effects: `intercept` (number) + `coefficients` (map contrast-column → β). |
| `random` | object | Random-effect structure by unit (`subject`, `item`). Empty `{}` ⇒ no random effects. |
| `response` | object | Outcome family + parameters (see below). |

### Factors

```json
{ "name": "condition",
  "levels": ["related", "unrelated"],
  "contrasts": { "cond": [-0.5, 0.5] },
  "vary_within": ["subject", "item"] }
```

* `contrasts` maps one or more **contrast-column names** to a numeric value per level
  (length = number of levels). Fixed coefficients and random slopes are keyed by these
  contrast-column names. This is exactly how effects are specified in `lme4`/DeBruine &
  Barr (2021): effects are coefficients on contrast-coded predictors.
* `vary_within`: the factor is crossed *within* the listed units (a within-unit factor),
  expanding each unit combination into one row per level.
* `between`: `"subject"` or `"item"` — the factor partitions that unit into equal blocks
  in level order (a between-unit factor; does not expand rows).

### Random effects (per unit)

```json
"subject": {
  "intercept_sd": 0.12,
  "slopes": { "cond": 0.04 },
  "correlations": { "intercept,cond": 0.2 }
}
```

The random-effect **column order** is `["intercept", <slopes in listed order>]`. A
covariance matrix `Σ = D · R · D` is formed from the SDs `D` and correlation matrix `R`
(unit diagonal; off-diagonals from `correlations`, keyed `"a,b"`). Per unit, a vector
`b = L z` is drawn, where `L` is the lower Cholesky factor of `Σ` and `z` are iid standard
normals. The unit's contribution to a row's linear predictor is
`b[intercept] + Σ_k b[slope_k] · (contrast value of slope_k for that row)`.

### Response families

| `family` | Parameters | Generation |
|---|---|---|
| `gaussian` | `sigma` | `y = η + σ·z` |
| `shifted_lognormal` | `sigma`, `shift` | `y = shift + exp(η + σ·z)` (reaction times) |
| `bernoulli` | — | `p = invlogit(η)`, `y = 1[u < p]` (accuracy; logit link) |
| `poisson` | — | `λ = exp(η)`, `y =` inverse-CDF Poisson (counts; log link) |
| `ordinal` | `thresholds` (K−1 cut-points) | cumulative-logit: `P(Y≤k) = invlogit(θ_k − η)` (Likert) |

`η` (the linear predictor for a row) = `intercept + Σ β_col · contrast_col +`
subject random part `+` item random part. `name` sets the output column name; optional
`round` sets decimal rounding of the response.

## RNG contract (identical across all implementations)

**Uniform generator.** L'Ecuyer (1988) combined LCG:

```
s1 ← (40014 · s1) mod 2147483563
s2 ← (40692 · s2) mod 2147483399
d  ← s1 − s2 ;  if d < 1 then d ← d + 2147483562
u  ← d / 2147483563            # u ∈ (0, 1)
```

All products stay below 2^53, so the arithmetic is exact in IEEE-754 doubles and in
Python integers alike. Seeding: `s1 ← 1 + (|seed| mod 2147483562)`,
`s2 ← 1 + ((40692 · s1) mod 2147483398)`, then 10 warm-up draws are discarded.

**Normal deviates.** Wichura (1988) Algorithm AS 241 applied to `u` (the algorithm R's
`qnorm` uses); deviates therefore agree to full double precision.

**Draw order (must be identical everywhere):**

1. For each subject `s = 1..S`: draw `q_subject` standard normals (intercept, then each
   slope in listed order); set `b_subject[s] = L_subject · z`.
2. For each item `t = 1..I` (if items exist): draw `q_item` standard normals; set
   `b_item[t] = L_item · z`.
3. Iterate observations in **canonical row order** and draw exactly one response deviate
   (normal for gaussian/lognormal; uniform for bernoulli/poisson/ordinal) per row.

**Canonical row order:** nested loops, outermost first —
`for s in 1..S: for t in 1..I: for (each within-factor level-combination, factors in
listed order, levels in listed order): emit row`. Between-unit factors assign a level to
each unit by equal blocks in level order and do **not** expand rows.

## References

* L'Ecuyer, P. (1988). Efficient and portable combined random number generators.
  *Communications of the ACM, 31*(6), 742–751.
* Wichura, M. J. (1988). Algorithm AS 241: The percentage points of the normal
  distribution. *Applied Statistics, 37*(3), 477–484.
* DeBruine, L. M., & Barr, D. J. (2021). Understanding mixed-effects models through data
  simulation. *Advances in Methods and Practices in Psychological Science, 4*(1).
