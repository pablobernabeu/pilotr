# The `pilotr` design specification (v0.3)

A design specification is a single JSON object describing a data-generating process (DGP) for
an experiment. It is the contract shared by the web application, the R package and the Python
package. Given the same spec and seed, every implementation must produce an identical data set,
within the boundary set out under [Scope of the guarantee](#scope-of-the-guarantee) below.

The machine-readable form of this document is
[`design.schema.json`](https://github.com/pablobernabeu/pilotr/blob/main/spec/design.schema.json).
The link is absolute because this file is read both here in the repository and embedded in the
documentation site, where a path relative to the repository would not resolve. Both
implementations enforce it through `validate_spec()`, together with the cross-field rules JSON
Schema cannot express.

## Versioning

`spec_version` is a `"major.minor"` string. A specification without it is a 0.2 specification,
which is what every specification written before the field existed is.

A specification that uses a feature introduced in 0.3 must declare `"0.3"` or later. This is not
bookkeeping: a 0.2 implementation reads such a specification differently and generates different
data while reporting success, which is the worst failure mode available to a reproducibility tool.
The 0.3 features are observation-level predictors, `dist`, `reliability`, the `exgaussian` family,
the `correlated` flag and interaction random slopes. An implementation from 0.3 onwards refuses a
specification declaring a version newer than it understands, and never reads it in part.

## Reading a specification

The same file must mean the same design to every implementation. JSON leaves a few things to the
reader, and R's jsonlite and Python's json settle them differently, so the rules below fix them.
Both implementations' `load_spec()` follow them, and `tools/parity/validate_cross.py` reads a
battery of awkward files through both to check that they agree.

A specification file is UTF-8, with or without a byte-order mark, whatever the locale of the
machine reading it, so level labels such as `"fácil"` or `"Łatwy"` read the same everywhere.
Every field is looked up by its exact name. R's `$` used to match a name partially, so a `random`
entry named `subject_site` was also read as `random.subject`.

A file that repeats a key within one object is refused. JSON leaves the meaning of a repeated key
undefined (RFC 8259, section 4), and the two readers differed: given `{"grp": 0.3, "grp": 0.5}` as
coefficients, R applied 0.3 twice and Python 0.5 once. The refusal happens as the file is read, so
it applies even when validation is switched off.

A JSON `null` in an optional field means the field is absent, so `"spec_version": null` is a 0.2
specification and `"slopes": null` means no slopes. A `null` in a required field counts as missing.
A unit given as `"item": null` is refused, since the object it stands for is not there, and so is a
`null` inside an array, such as a level.

Where the schema has a single value, a one-element array holding that value is read as the value,
so `"seed": [2024]` is the seed 2024. That is how `jsonlite::write_json()` writes every scalar by
default, and such a file reads in Python as it does in R. For the same reason, `[]` is read as an
empty object at exactly `random`, `fixed.coefficients`, `random.<group>.slopes` and
`random.<group>.correlations`. Anywhere else an object belongs, such as `units`, `fixed`,
`response` or `contrasts`, `[]` is refused. In the other direction, a `vary_within` given as one
string and a single threshold given as a bare number are read as one-element arrays, because
pilotr's own `spec_json()` wrote both forms before 0.3. A count written with a decimal point, such
as `"n": 8.0`, is the integer 8, as JSON Schema draft-07 counts it.

jsonlite turns an array of mixed types, such as levels `["a", 1]`, into strings, and R then
accepts it. That reading is R's alone: Python refuses such an array, and writing every level as a
string avoids the difference.

A seed must be a whole number between −9007199254740991 and 9007199254740991, that is
±(2^53 − 1). jsonlite reads a JSON number as a double, which holds every integer in that range and
not every one beyond it, while Python reads an integer of any size exactly. A larger seed was
therefore two different numbers in the two readers, which drew different data.

## Top-level fields

A specification is a flat object with the nine fields below, of which `name`, `seed`, `units`,
`fixed` and `response` are required. The rest may be omitted: a specification with no
`spec_version` is a 0.2 specification, and one with no `factors`, `predictors` or `random` simply
has none of them.

| Field | Type | Meaning |
|---|---|---|
| `spec_version` | string | Specification version, `"major.minor"`. Absent means 0.2. |
| `name` | string | Human label for the design. |
| `seed` | integer | Master seed, within ±(2^53 − 1) (see RNG contract below). |
| `units` | object | Sampling units, e.g. `{"subject": {"n": 30}, "item": {"n": 24}}`. `item` is optional. Add `per_subject` to `item` (e.g. `{"n": 40, "per_subject": 12}`) for partial crossing, in which each subject sees a random subset of items. |
| `factors` | array | Experimental factors (categorical; see below). |
| `predictors` | array | Optional continuous predictors (see below). |
| `fixed` | object | Fixed effects: `intercept` + `coefficients` (map column → β). A coefficient key may be a single column or an `"a:b"` interaction (the product of columns a and b). |
| `random` | object | Random-effect structure by unit (`subject`, `item`). Empty `{}` ⇒ no random effects. |
| `response` | object | Outcome family + parameters (see below). |

A coefficient or slope key that names no existing column contributes zero, so it silently removes
the term and raises nothing. `validate_spec()` therefore refuses such a key: a design whose focal
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
* `vary_within`: the factor varies *within* units (a within-unit factor), expanding each unit
  combination into one row per level. pilotr crosses a within factor with every unit of the
  design, whatever the list names, so in a design with items the list is `["subject", "item"]`.
  A list that leaves a unit out has never changed the data. It is deprecated and draws a
  warning, and from spec version 0.4 it will be refused. A factor whose levels are carried by
  the items, each item appearing at one level only, is `"between": "item"` instead. Since a
  within factor is crossed with every unit, each subject–item pair of the design is observed at
  every level. A counterbalanced design, which shows each subject each item at one level only,
  is written with a list factor between subjects and an item-set factor between items, as in the
  [worked encodings](#worked-encodings) below.
* `between`: `"subject"` or `"item"`. The factor assigns each unit one level and does not
  expand rows. With `N` units and `L` levels, unit `u` (counted from 1) takes the level with
  index `⌊(u − 1)·L/N⌋` (counted from 0), so the units fall into `L` consecutive blocks in level
  order. The blocks are equal when `N` is a multiple of `L` and otherwise differ in size by one.

A factor sets exactly one of `vary_within` and `between`, and lists each level once.

Each between factor is allocated by that rule on its own, so two between factors over the same
unit would fall into the same or overlapping blocks. Over 40 subjects, two two-level factors
would give cells of 20, 0, 0 and 20, and neither the second effect nor the interaction could be
estimated. A specification in which two or more factors are between the same unit is therefore
refused. A factorial between design is written as one factor whose levels are the cells, as in
the first of the [worked encodings](#worked-encodings). Specification version 0.4 will allocate
several between factors jointly.

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
was silently given one value per item. It is now validated against the three names that exist.

`dist` (new in 0.3) selects the distribution, either `"normal"`, which uses `mean` and `sd`, or
`"uniform"`, which uses `min` and `max`. A uniform draw consumes exactly as much of the random
stream as a normal one, since a normal is produced by transforming a single uniform, so switching
between them does not move the stream.

#### Reliability

`reliability` (new in 0.3) simulates a predictor measured with error:

```json
{ "name": "z_reading", "varies_by": "subject", "sd": 1, "reliability": 0.8 }
```

The latent value drives the linear predictor and any random slope keyed on the predictor, while
the returned data carries the observed, contaminated value, which is what an analyst would have
measured. Writing `ρ` for the reliability and using the predictor's population
mean and standard deviation,

```
observed = mean + (true − mean + sd·sqrt((1 − ρ)/ρ)·z) · sqrt(ρ)
```

so the observed variable has the same variance as the latent one and correlates `sqrt(ρ)` with it.
Reliability in the classical sense is that squared correlation, which is why the field is `ρ`.

The attenuation is `sqrt(ρ)`, where the textbook regression-dilution result gives `ρ`, because
both variables are placed on the same variance here. Standardising the observed variable back to
the latent one's variance absorbs the `1/sqrt(ρ)` factor that result carries.

The moments used are the population ones. R's `mean()` and `sd()` accumulate in long double and
Python's do not, so standardising against the sample mean and standard deviation of the values
drawn would reintroduce a cross-language divergence.

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
separator, `"a~b"`, is also accepted). Each element is computed as `(sd_i · sd_j) · r_ij`. That
bracketing is part of the contract, because floating-point multiplication is not associative and
the alternative grouping lands on a different double for some inputs. Per unit, a vector
`b = L z` is drawn, where `L` is the lower Cholesky factor of `Σ` and `z` are iid standard
normals. The unit's contribution to a row's linear predictor is
`b[intercept] + Σ_k b[slope_k] · (design value of slope_k for that row)`.

Slope keys follow exactly the same rule as fixed coefficients: a contrast column, a continuous
predictor or an `"a:b"` interaction between them. Before 0.3 an interaction slope was accepted,
sized into the covariance, drawn from the stream and then silently discarded, so the emitted
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
The units are assigned to groups by the block rule of a between factor. With `N` units and `K`
groups, unit `u` joins group `⌊(u − 1)·K/N⌋`, counted from 0 and written to the data counted
from 1. The groups are therefore equal when `N` is a multiple of `K`. For example, subjects
nested in clusters:

```json
"site": { "over": "subject", "n": 12, "intercept_sd": 0.5, "slopes": { ... } }
```

Each group draws a random-effect vector (intercept + any slopes) applied to all rows of the
units in that group, and the simulated data gains a column with the group id. Useful for
hierarchical designs (e.g. participants within sites, schools or languages).

A between factor and a grouping factor over the same unit follow the same rule, so they nest.
When `K` is a multiple of the factor's `L` levels, each group lies wholly in one level. With 120
subjects in 12 sites and a two-level between factor, sites 1 to 6 are all in the first
condition and sites 7 to 12 in the second, which randomises whole clusters. When `L` is a
multiple of `K`, each level lies wholly in one group, and otherwise a group can straddle two
levels. A design randomised within clusters, each cluster holding both conditions, is written
with a between factor whose levels nest in the groups, as in the last of the
[worked encodings](#worked-encodings).

### Worked encodings

Three common designs have no field of their own and are written with the fields above. Each is
also a case of the parity harness in `tools/parity/cases/`, so its data are held bit-identical
across the two implementations and anchored by a recorded hash.

A 2 × 2 between-subjects design is one between factor whose four levels are the cells. Its
contrast columns carry the two main effects, and the coefficient keyed `a:b` the interaction, as
in [`between_cells_2x2.json`](https://github.com/pablobernabeu/pilotr/blob/main/tools/parity/cases/between_cells_2x2.json):

```json
"factors": [
  { "name": "cell",
    "levels": ["a1.b1", "a1.b2", "a2.b1", "a2.b2"],
    "contrasts": { "a": [-0.5, -0.5, 0.5, 0.5], "b": [-0.5, 0.5, -0.5, 0.5] },
    "between": "subject" }
],
"fixed": { "intercept": 10, "coefficients": { "a": 0.5, "b": 0.3, "a:b": 0.2 } }
```

The file's 40 subjects fall into four blocks of 10, one per cell. A factor with more levels
takes one contrast column per degree of freedom, and its interaction with another factor one key
per pair of columns. `spec_from_model()` in R writes this encoding itself when it reads two
factors between the same unit off a fitted pilot.

A two-list counterbalanced design shows each subject each item once, at one of two levels, and
swaps the levels between two lists of subjects. It is written with a list factor between
subjects and an item-set factor between items, whose contrasts multiply to the condition. With
`l` at −1 and 1 and `g` at −0.5 and 0.5, the product `l·g` is −0.5 or 0.5, so the condition
effect and its random slopes are keyed `l:g`, as in
[`two_list_counterbalanced.json`](https://github.com/pablobernabeu/pilotr/blob/main/tools/parity/cases/two_list_counterbalanced.json):

```json
"factors": [
  { "name": "list", "levels": ["list1", "list2"], "contrasts": { "l": [-1, 1] },
    "between": "subject" },
  { "name": "item_set", "levels": ["setA", "setB"], "contrasts": { "g": [-0.5, 0.5] },
    "between": "item" }
],
"fixed": { "intercept": 6, "coefficients": { "l:g": 0.06 } },
"random": {
  "subject": { "intercept_sd": 0.12, "slopes": { "l:g": 0.04 } },
  "item": { "intercept_sd": 0.08, "slopes": { "l:g": 0.02 } }
}
```

Random slopes keyed on an interaction are a 0.3 feature, so the file declares
`"spec_version": "0.3"`. With 24 subjects and 18 items it gives 432 rows, one for each
subject–item pair, half the 864 of the fully crossed design with a within factor. Every subject
sees nine items in each condition, and every item appears in both conditions across the two
lists. In R, `model_formula()` gives `.y ~ l_g + (1 + l_g | subject) + (1 + l_g | item)`, where
`l_g` is the condition column that `model_data()` builds. The encoding suits two conditions,
since with three or more each condition is a combination of several product columns and the
analysis model has to be written by hand.

Randomisation within clusters, with both conditions in every cluster, is written with a between
factor whose levels nest in the clusters. With 120 subjects in 12 sites, a factor of 24 levels
gives each level a block of five subjects and each site two levels. A contrast that alternates
across the levels then puts five subjects in each condition at every site, as in
[`within_cluster_randomised.json`](https://github.com/pablobernabeu/pilotr/blob/main/tools/parity/cases/within_cluster_randomised.json):

```json
"factors": [
  { "name": "arm",
    "levels": ["s01_control", "s01_treatment", "s02_control", ..., "s12_treatment"],
    "contrasts": { "trt": [-0.5, 0.5, -0.5, ..., 0.5] },
    "between": "subject" }
],
"random": {
  "site": { "over": "subject", "n": 12, "intercept_sd": 3, "slopes": { "trt": 1 } }
}
```

The labels name the site each level falls in, which holds because 24 is a multiple of 12. Since
every site holds both conditions, the by-site slope on `trt`, a treatment effect that varies
between sites, can be estimated.

### Response families

`response.family` selects one of eight generation rules, each mapping the linear predictor `η` to
an outcome on the scale that family works on. The parameters column names the extra fields the
family reads from `response`.

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

The Poisson mean has an upper bound, and it is the same in every implementation. The inverse-CDF
walk starts from `exp(−λ)`, which underflows to exactly zero once `λ` passes about 746, and from
there the cumulative distribution can never reach the drawn uniform, so no count exists to return.
Both implementations refuse such a mean with the same message, where earlier versions returned the
walk's iteration cap as though it were a draw. The effect is to cap the linear predictor of a
`poisson` response at roughly 6.6, a mean count near 750, which is well above the rates count
outcomes are normally specified at.

`η` (the linear predictor for a row) = `intercept + Σ β_key · value(key)`, where a key is a
contrast column, a continuous predictor or an `"a:b"` interaction (the product of the named
columns), `+` subject random part `+` item random part `+` the random parts of any additional
grouping factors. `name` sets the output column name. An optional `round` sets the decimal
rounding of the response, and applies only to the families whose outcome is continuous, the others
being integers already.

### Names

The names in a specification become columns, of the simulated data or of the analysis data that
`model_data()` builds from it. `validate_spec()` refuses a specification in which two of them
would land in one place. Each case below used to validate and then move an effect, rescale a
variance or overwrite a column without a word. R replaced a column written twice where Python
appended a second one under the same name, so the two implementations also exported different
tables from one specification.

1. The simulated data has these columns, in this order: `subject`, `item` when the design has
   items, one per additional grouping factor, one per factor, one per continuous predictor and
   the response. Their names are non-blank and all different.
2. A contrast column belongs to one factor. Its name differs from every predictor, from
   `subject`, from `item` when the design has items, from every additional grouping factor, from
   every other factor and from the response. It may carry its own factor's name, as two-level
   designs often do.
3. The analysis column `a_b` that `model_data()` builds for an interaction key `a:b` does not
   take the name of another column, another interaction's included.
4. A factor lists each level once.
5. A factor sets exactly one of `vary_within` and `between`.
6. A correlation key names two different terms, and each pair is given once, in one order.
7. `random.item` appears only in a design with an item unit. Without one it was dropped from the
   simulated data, while `model_formula()` still put an item term into the analysis.

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
`**`. R special-cases small integer exponents while Python calls the library `pow()`, and
measured over 200,000 draws in the Gamma sampler's range the two disagreed on a third of inputs
by up to 6 ulp. Because that value decides a rejection step, the disagreement also changed how
many draws were consumed.

## RNG contract (identical across all implementations)

Two generators and one draw order fix the random stream. Everything drawn anywhere in pilotr
comes from the uniform generator below, either directly or through the normal transform, and in
the sequence set out under [Draw order](#draw-order).

### Uniform generator

Uniform deviates come from L'Ecuyer's (1988) combined LCG:

```
s1 ← (40014 · s1) mod 2147483563
s2 ← (40692 · s2) mod 2147483399
d  ← s1 − s2 ;  if d < 1 then d ← d + 2147483562
u  ← d / 2147483563            # u ∈ (0, 1)
```

All products stay below 2^53, so the arithmetic is exact in IEEE-754 doubles and in
Python integers alike. The seeding rule is `s1 ← 1 + (|seed| mod 2147483562)` and
`s2 ← 1 + ((40692 · s1) mod 2147483398)`, after which 10 warm-up draws are discarded. The
seed is bounded by ±(2^53 − 1), so it too is held exactly in both, and the remainder is taken of
the same integer (see [Reading a specification](#reading-a-specification)).

### Normal deviates

Every normal in pilotr is Wichura's (1988) Algorithm AS 241 applied to `u`, the algorithm R's
`qnorm` uses, so the two languages agree to full double precision.

### Draw order

The order in which the stream is consumed is as much part of the contract as the generator
itself, and every implementation must follow the sequence below exactly.

0. If `units.item.per_subject` is set (partial crossing): for each subject `s = 1..S`, in
   row-build order, sample that subject's item subset by a partial Fisher–Yates shuffle,
   consuming one uniform per sampled item (`per_subject` uniforms per subject). These are
   the first RNG draws. (Skipped entirely under full crossing, so fully crossed specs keep
   the original stream.)
1. For each continuous predictor (in listed order): for each of its units `u = 1..N`, draw one
   deviate, `N(mean, sd)` by default or `Uniform(min, max)` when `dist` is `"uniform"`, and then,
   only when `reliability` is present and below 1, one further standard normal for that value's
   measurement error. `N` is the number of subjects, of items or of rows, according to
   `varies_by`. (Skipped entirely when there is no `predictors` block, so factor-only specs keep
   the original stream. A uniform costs the same one draw as a normal, and a `reliability` of 1 or
   absent costs nothing, so neither moves the stream either.)
2. For each subject `s = 1..S`: draw `q_subject` standard normals (intercept, then each
   slope in listed order), and set `b_subject[s] = L_subject · z`.
3. For each item `t = 1..I` (if items exist): draw `q_item` standard normals, and set
   `b_item[t] = L_item · z`.
4. For each additional grouping factor (in the order the `random` entries are listed): for
   each group `g = 0..K−1`, draw `q_group` standard normals, and set `b_group[g] = L_group · z`.
5. Iterate observations in canonical row order, defined below, and draw the response. Most
   families consume exactly one deviate per row, a normal for gaussian, lognormal and
   shifted_lognormal and a uniform for bernoulli, poisson and ordinal. The ex-Gaussian consumes
   exactly two, a normal then a uniform, in that order. The beta family instead consumes a
   variable, data-dependent number of draws per row: two Gamma variates through the
   Marsaglia–Tsang rejection sampler, each consuming normal–uniform pairs until acceptance (with
   one extra uniform per Gamma variate whose shape is below 1).

### Canonical row order

Step 5 iterates the observations as nested loops, outermost first,
`for s in 1..S: for t in 1..I: for (each within-factor level-combination, factors in
listed order, levels in listed order): emit row`. A between factor assigns each unit a level by
the block rule under [Factors](#factors) and does not expand rows.

### Extending the draw order

A new feature must consume no draws at all when it is not used, so that every specification
written before it stays bit-identical. Each of the 0.3 additions is built that way: an
observation-level predictor only appears when declared, a uniform draw costs the same as the
normal it replaces, a `reliability` of 1 or absent draws nothing, and a family branch is
reached only when that family is selected.

### Replicate seeds

The power and precision loops derive their per-replicate seeds from the specification's seed by
drawing from the shared generator and skipping any duplicate. Adding the replicate index, the
rule they used before 0.3, does not give independent streams here: seeding sets `s1` to
`1 + (seed mod 2147483562)` and `s2` from `s1`, with only ten warm-up draws discarded, so the
first draw of replicate `i` correlated 0.95 with that of replicate `i + 1`. An arithmetic
scramble does not help, since the seeding rule is linear in the seed. This changed every number
pilotr produced before 0.3.

## Scope of the guarantee

Identical data means bit-identical, and holds exactly for the `gaussian` family and for any design
whose response path applies no transcendental function to the linear predictor. It also holds for
every family when `response.round` is set, since rounding quantises away a last-bit difference.

For `lognormal`, `shifted_lognormal`, `exgaussian`, `bernoulli`, `poisson`, `ordinal` and `beta`
without `round`, results may differ in the last unit in the last place. IEEE-754 requires correct
rounding for addition, subtraction, multiplication, division and square root, but leaves the
rounding of `exp()` and `log()` to the implementation. The R and Python builds on a given
platform need not share a maths library. Measured over 200,000 arguments in the log-reaction-time
range, R and CPython `exp()` disagreed on 0.44% of them by up to 6 ulp, and `log()` on 0.12% by up
to 1 ulp.

This is demonstrable, and has been demonstrated. Taking a shifted-lognormal design and switching
only its family to `gaussian`, so that the seed, the random-effect structure, the linear predictor
and the entire draw sequence are unchanged, gives bit-identical output in both languages, while the
lognormal original differs in a handful of rows. `exp()` is the only remaining difference.

For the discrete families the practical consequence is different in kind. A last-bit difference in
`exp()` usually changes nothing, because the outcome is an integer decided by a comparison. When
the comparison sits exactly on a threshold, though, the outcome moves by a whole category. That is
rare but possible, so a design analysis that has to be reproducible to the last observation should
either set `round` or stay with `gaussian`.

The stricter guarantee within one language is unconditional: the same implementation, specification
and seed always produce the same data.

pilotr simulates no attrition, exclusion or missing response, so every data set is complete and
the sample size in a specification is the number of units analysed. A study that expects to lose a
proportion `p` of its participants recruits `N/(1 − p)` to analyse `N`, so 60 analysed at 10%
attrition means recruiting 67. Missingness that depends only on what is observed, such as the
condition a participant was assigned to, is ignorable for likelihood-based inference (Rubin,
1976), so it costs power and balance without biasing the estimates. Missingness that depends on
the values that would have been observed, which pilotr does not simulate, can bias them as well.

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
* Rubin, D. B. (1976). Inference and missing data. *Biometrika, 63*(3), 581–592.
  https://doi.org/10.1093/biomet/63.3.581
