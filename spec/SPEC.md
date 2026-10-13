# The `pilotr` design specification (v0.4)

A design specification is a single JSON object describing a data-generating process (DGP) for
an experiment. It is the contract shared by the web application, the R package and the Python
package. Given the same spec and seed, every implementation must produce an identical data set,
within the boundary set out under [Scope of the guarantee](#scope-of-the-guarantee) below.

The machine-readable form of this document is
[`design.schema.json`](https://github.com/pablobernabeu/pilotr/blob/main/spec/design.schema.json).
The link is absolute because this file is read both here in the repository and embedded in the
documentation site, where a path relative to the repository would not resolve. The schema states
which fields each object may hold and their types, ranges and allowed values. It also states
which parameters each response family requires and refuses. Both implementations enforce those
rules through `validate_spec()` in strict mode, its default, and `tools/parity/validate_cross.py`
checks that the schema accepts every specification of its battery that both validators accept.
Of the readings set out under [Reading a specification](#reading-a-specification), the schema
accepts the single-value forms of `vary_within` and `thresholds` and a count written as `8.0`. It
refuses a `null` standing for an absent field, `[]` standing for an empty object and a
one-element array standing for a single value, all of which both implementations accept.

`validate_spec()` also enforces the rules that relate one value to another or to the
implementation reading the file. Every coefficient and slope key names a contrast column or a
predictor, a contrast has one value per level and thresholds increase. A uniform predictor's
`min` lies below its `max`, `per_subject` does not exceed the number of items, and a unit that a
factor, a predictor, a grouping factor or `random.item` names exists. A correlation key names two
terms of its group, and each pair is given once. A specification that uses a feature introduced
in 0.3 or 0.4 declares that version or later, and none may declare a version newer than the
implementation. Names do not collide (see [Names](#names)). Below version 0.4, at most one factor
is between each unit, and from 0.4 a within factor lists every unit of the design (see
[Factors](#factors)). Lenient validation (`strict = FALSE`) reports an unknown field as a warning
where the schema refuses it.

## Versioning

`spec_version` is a `"major.minor"` string. A specification without it is a 0.2 specification,
which is what every specification written before the field existed is.

A specification that uses a feature introduced in 0.3 must declare `"0.3"` or later. This is not
bookkeeping: a 0.2 implementation reads such a specification differently and generates different
data while reporting success, which is the worst failure mode available to a reproducibility tool.
The 0.3 features are observation-level predictors, `dist`, `reliability`, the `exgaussian` family,
the `correlated` flag and interaction random slopes. An implementation from 0.3 onwards refuses a
specification declaring a version newer than it understands, and never reads it in part.

Version 0.4 adds one feature, several between factors on one unit, which it allocates jointly (see
[Factors](#factors)). Up to 0.3, each between factor was allocated on its own, so two over one
unit fell into the same or overlapping blocks. A specification that has them and declares an
earlier version is therefore refused. Version 0.4 also refuses a within factor whose
`vary_within` leaves out a unit of the design, which draws only a warning in a specification
declaring an earlier version.

An implementation that writes a specification declares the lowest version the specification
needs, and never one below 0.3. The R package's `build_spec()` and `spec_from_model()` therefore
declare 0.3 unless the specification uses a 0.4 feature, so the releases that read 0.3, R 0.3.1
and Python 0.3.0 among them, read what they write. Declaring 0.2 would gain no reader, since
releases before 0.3 ignore the field.

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

Both readers convert a decimal with correct rounding, jsonlite through the C library's `strtod()`
and Python through `float()`, so a number means the same double to both. A writer keeps a value
exactly by giving enough significant digits for correct rounding to return it, and 17 always
suffice. R's `spec_json()` writes the shortest such form, and `tools/parity/json_roundtrip.py`
checks that Python reads back the doubles R wrote.

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
  A list that leaves a unit out has never changed the data. It draws a warning in a
  specification declaring a version below 0.4 and is refused from 0.4. A factor whose levels
  are carried by the items, each item appearing at one level only, is `"between": "item"`
  instead. Since a within factor is crossed with every unit, each subject–item pair of the
  design is observed at every level. A counterbalanced design, which shows each subject each
  item at one level only, is written with a list factor between subjects and an item-set factor
  between items, as in the [worked encodings](#worked-encodings) below.
* `between`: `"subject"` or `"item"`. The factor assigns each unit one level and does not
  expand rows. With `N` units and `L` levels, unit `u` (counted from 1) takes the level with
  index `⌊(u − 1)·L/N⌋` (counted from 0), so the units fall into `L` consecutive blocks in level
  order. The blocks are equal when `N` is a multiple of `L` and otherwise differ in size by one.
  Several between factors on one unit are allocated together, as set out below.

A factor sets exactly one of `vary_within` and `between`, and lists each level once.

From version 0.4, the between factors of one unit are allocated jointly. Taken in listed order,
with `L_1, …, L_k` levels, they have `C = L_1·…·L_k` combinations of levels, the cells. Unit `u`
of `N` takes cell `⌊(u − 1)·C/N⌋`, counted from 0, so the units fall into `C` consecutive blocks.
The blocks are equal when `N` is a multiple of `C` and otherwise differ in size by one. A cell is
read as one level per factor, with the first-listed factor varying slowest and the last fastest.
With two factors, the second of three levels, cell `c` holds level `⌊c/3⌋` of the first factor and
level `c mod 3` of the second. With one factor, this is the rule above. Over 40 subjects, two
two-level factors give four cells of 10. The allocation is integer arithmetic and draws nothing.

Up to version 0.3, each between factor was allocated by the rule above on its own, so two between
factors over the same unit fell into the same or overlapping blocks. Over 40 subjects, two
two-level factors gave cells of 20, 0, 0 and 20, and neither the second effect nor the interaction
could be estimated. Two or more factors between the same unit are therefore refused in a
specification declaring a version earlier than 0.4. A factorial between design is written either
as several between factors in a 0.4 specification or as one factor whose levels are the cells,
which every version reads, as in the first of the [worked encodings](#worked-encodings). The two
allocate the units alike and, given the same seed, draw the same responses.

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
levels. The cells of several between factors on the unit nest in the same way, their number `C`
taking the part of `L`. With 12 sites and a 2 × 2 design, each site lies wholly in one of the four
cells. A design randomised within clusters, each cluster holding both conditions, is written
with a between factor whose levels nest in the groups, as in the last of the
[worked encodings](#worked-encodings).

### Worked encodings

Three common designs are written with the fields above in forms that need no version above 0.3.
Each is also a case of the parity harness in `tools/parity/cases/`, so its data are held
bit-identical across the two implementations and anchored by a recorded hash.

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

From version 0.4, the same design can also be written as two between factors, `A` with contrast
`a` and `B` with contrast `b`, as in
[`two_between_joint.json`](https://github.com/pablobernabeu/pilotr/blob/main/tools/parity/cases/two_between_joint.json).
Its subjects fall into the same four cells in the same order, and with the same seed it draws the
same scores. The two data sets differ only in their label columns, which are `A` and `B` in one
and `cell` in the other. A release that reads no version above 0.3 reads only the encoding with
one factor, which is why `spec_from_model()` keeps it.
[`three_by_two_joint.json`](https://github.com/pablobernabeu/pilotr/blob/main/tools/parity/cases/three_by_two_joint.json)
crosses a 3 × 2 design between 12 subjects, two to a cell, with a 2 × 2 design between 8 items
and a within factor.

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
walk starts from `exp(−λ)`, which is subnormal once `λ` exceeds `−log` of the smallest normal
double, `1022·ln 2`, about 708.4. A subnormal carries fewer significant bits. From a mean of
about 738, the walk's cumulative distribution no longer sums to 1, so its counts are biased or the
drawn uniform lies out of its reach. Beyond about 745.1, `exp(−λ)` underflows to zero. Every
implementation therefore refuses a mean above 708.3964185322641, the double nearest `1022·ln 2`,
before the walk starts, and with the same message. Earlier versions served means up to about 745,
with biased counts from about 738. The effect is to cap the linear predictor of a `poisson`
response at about 6.56, a mean count of about 708, which is well above the rates count outcomes
are normally specified at. The walk itself is set out under
[Accumulation order](#accumulation-order-identical-across-all-implementations).

A `lognormal` or `shifted_lognormal` response is refused when `exp(η + σ·z)` overflows, as it does
once its argument passes about 709.78. The message is the same in every implementation and names
`η`, since `exp(η)` can be finite in a row whose normal draw carries the argument past the limit.
Earlier versions wrote an infinite response in R and raised a bare overflow error in Python.

A `beta` response lies strictly inside (0, 1), the support of the Beta distribution and the range a
Beta regression accepts. It is `x / (x + y)`, where `x` is a Gamma variate with shape `μ·φ` and
`y`, drawn after it, one with shape `(1−μ)·φ`. Two cases are refused before the division, with the
same messages in every implementation. The first is a shape of exactly zero. It arises when
`invlogit(η)` rounds to 1 (from an `η` of about 36.74) or underflows to 0 (below about −745.1), or
when `φ` is small enough to take either product to zero by itself. The second is a pair of
variates that both underflow to zero, leaving the ratio undefined, as can happen when both shapes
are tiny. Earlier versions returned 0, 1 or `NaN` in R and raised a division error in Python.

Otherwise, the ratio rounds to exactly 1 whenever `y` is below about 2^−54 of `x`. That happens at
realistic means near a ceiling, as in 161 of 4,000 rows of a design with an intercept of 3 and
`φ = 2`. The ratio is exactly 0 when `x` underflows. A ratio of exactly 1 becomes 1 − 2^−53, the
largest double below 1, and one of exactly 0 becomes 2^−1074, the smallest positive double. Where
both variates were representable, a ratio that rounded to 1 lay between 1 − 2^−53 and 1, so no
double inside the interval is nearer the draw. The replacement is a comparison and draws nothing,
so the stream is unchanged. When `μ·φ` or `(1−μ)·φ` is well below 1, part of the distribution lies
nearer a boundary than doubles resolve, and those rows sit at the two bounds. A rounded `beta`
response is kept inside the interval too, by the last rule under
[Rounding](#rounding-identical-across-all-implementations).

`η` (the linear predictor for a row) = `intercept + Σ β_key · value(key)`, where a key is a
contrast column, a continuous predictor or an `"a:b"` interaction (the product of the named
columns), `+` subject random part `+` item random part `+` the random parts of any additional
grouping factors. `name` sets the output column name. An optional `round` sets the decimal
rounding of the response, and applies only to the families whose outcome is continuous, the others
being integers already. For `beta` it must be at least 1, since every proportion rounded to a whole
number is 0 or 1. The rounding follows the rule under
[Rounding](#rounding-identical-across-all-implementations).

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

The Poisson walk is fixed in the same way. It starts from `p_0 = exp(−λ)` and `c_0 = p_0`, and
step `k` computes `p_k = (p_{k−1}·λ)/k`, multiplying before it divides, and then
`c_k = c_{k−1} + p_k`. The draw is the first `k` at which `u ≤ c_k`. Bracketed as
`p_{k−1}·(λ/k)`, most of the cumulative sums land on a different double, and a uniform falling
between the two versions of a sum gives a different count. At `λ = 3.7`, `c_3` is
0.4941532441504185 with the first bracketing and 0.49415324415041845 with the second. A uniform
of 0.4941532441504185 therefore draws 3 with the first and 4 with the second.

## Rounding (identical across all implementations)

Every implementation applies `response.round` by one algorithm, that of R's `round()` from R 4.0.0
onwards (`src/nmath/fround.c`). The R package calls `round()` itself, and the Python package ports
it as `pilotr.core.fround()`. Every arithmetic operation below is carried out in double precision
and rounded to nearest, as IEEE-754 requires. For a response `y` and `round = d`:

1. If `y` is zero, infinite or not a number, it is returned unchanged.
2. If `d` is 0, `y` becomes the nearest whole number, a tie going to the even one.
3. Otherwise, let `x = |y|`, with binary exponent `e`, so that `2^e ≤ x < 2^(e+1)`. If
   `log10(2) · (e + 0.5) + d > 15`, the rounding position lies beyond the 15 significant digits a
   double always carries, and `y` is returned unchanged.
4. Otherwise, let `p = 10^d`, `t = x · p`, `lo = floor(t) / p` and `hi = ceil(t) / p`.
5. The result is `hi` if `hi − x < x − lo`, or if the two are equal and `floor(t)` is odd, and
   `lo` otherwise, with the sign of `y`.

`10^d` is computed by binary powering, as R's `R_pow_di()` computes it. It works through the
binary digits of `d` from the lowest, multiplying the result by the base wherever the digit is 1
and squaring the base before each further digit. The method matters from `d = 25`, where binary
powering and repeated multiplication land on different doubles, while up to `d = 22` the power is
exact either way. Beyond 308 places, where `10^d` would overflow, `x` is multiplied by `10^308` and
then by `10^(d − 308)`, and `floor(t)` and `ceil(t)` are divided by the two in the same order.
Beyond 323 places, `y` is returned unchanged. Validation accepts only a whole `d` of 0 or more.
Without it, a `d` that is not whole counts as `floor(d + 0.5)` from step 3 on, a negative `d`
takes `p = 1 / 10^|d|` and a `d` below −308 gives 0.

A `beta` response takes one further step, which keeps the result inside (0, 1). Let `q` be 1
multiplied by 10, `d` times, which is exact up to `d = 22`, and let `m = 1/q` for a `d` from 1 to
22 and 2^−1074 otherwise. A result of 0 or less becomes `m`. A result of 1 or more becomes `1 − m`,
or 1 − 2^−53 where `1 − m` rounds to 1, as it does from `d = 17`. Up to 22 places, `m` is the
double that step 5 returns for 10^−d, and up to 15 places, `1 − m` is the one it returns for
1 − 10^−d. A proportion rounded onto 0 or 1 therefore moves to the nearest value inside the
interval at `d` decimals, 0.001 or 0.999 at three. Validation refuses `d = 0` for `beta`, where
no such value exists. Before this rule, rounding sent 17% of the rows of one realistic low-mean
design to exactly 0. The parity case `tools/parity/cases/beta_boundary_rounded.json` checks the
rule end to end.

CPython's `round()` follows a different rule. It rounds the exact binary value of `y`, whereas
step 5 compares the two candidates by their distances from `x` as computed in double precision.
Within a few units in the last place of a decimal tie, the two rules can therefore take different
candidates. `round(1.05, 1)` is 1 in R, because in double precision `1.1 − 1.05` and `1.05 − 1`
are equal and `floor(t)`, 10, is even. It is 1.1 in CPython, because the double nearest 1.05 lies
just above it. `tools/parity/round_probes.R` writes 5,000 such values with R's results to
`python/tests/fixtures/r_round_probes.csv`, and CPython's `round()` took the other candidate for
166 of them. It agreed with R on all 1,000 random values in the same file. The parity case
`tools/parity/cases/round_near_tie.json` checks the rule end to end.

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

Every normal in pilotr is Wichura's (1988) Algorithm AS 241, the routine PPND16, applied to `u`.
It is the algorithm R's `qnorm` uses, and pilotr's R port returned `qnorm`'s double for each of
1,000,000 uniforms measured. With `q = u − 0.5` it has two branches:

```
if |q| ≤ 0.425:  r ← 0.180625 − q·q
                 z ← (q·A(r)) / B(r)
otherwise:       r ← u if q < 0, else 1 − u
                 r ← sqrt(−log(r))
                 z ← C(r − 1.6) / D(r − 1.6) if r ≤ 5, else E(r − 5) / F(r − 5)
                 z ← −z if q < 0
```

`A` to `F` are the paper's polynomials of degree 7, each evaluated in Horner form from its
highest power, `(…((c₇·r + c₆)·r + c₅)…)·r + c₀`, with the coefficients as `core.R` and `core.py`
write them. The central branch is arithmetic alone. The tail branch, which about 15% of draws
take, calls `log()`, whose last bit the maths library decides. Two implementations that do not
share that library can therefore turn the same uniform into different normals (see
[Scope of the guarantee](#scope-of-the-guarantee)).

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
   Marsaglia–Tsang rejection sampler. Each proposal draws a normal and, only when that normal
   gives `v > 0`, a uniform for the acceptance test. A variate whose shape is below 1 draws one
   further uniform, after the variate it boosts. The sampler is set out under
   [Sampling algorithms](#sampling-algorithms-identical-across-all-implementations).

### Canonical row order

Step 5 iterates the observations as nested loops, outermost first,
`for s in 1..S: for t in 1..I: for (each within-factor level-combination, factors in
listed order, levels in listed order): emit row`. The level combinations run with the last
within factor varying fastest. Under partial crossing, a subject's items run in ascending order.
A between factor assigns each unit a level by the block rule under [Factors](#factors), jointly
with any other between factor on the same unit, and does not expand rows.

### Extending the draw order

A new feature must consume no draws at all when it is not used, so that every specification
written before it stays bit-identical. Each of the 0.3 additions is built that way: an
observation-level predictor only appears when declared, a uniform draw costs the same as the
normal it replaces, a `reliability` of 1 or absent draws nothing, and a family branch is
reached only when that family is selected. The joint allocation of 0.4 draws nothing either, and
with one between factor per unit it allocates as earlier versions did.

### Replicate seeds

The power and precision loops derive their per-replicate seeds from the specification's seed. A
generator seeded with it draws `u`, the next seed is `⌊u · 2147483562⌋ + 1`, and a seed already
drawn is skipped, so no two replicates share one.

Adding the replicate index, the rule the loops used before 0.3, does not give independent
streams. Seeding is affine in the seed and each component recurrence is linear, so the `k`-th
draws of consecutive seeds differ by very nearly the same amount mod 1 (to within 1e−7). For the
first draw, the step `d` is −0.0091. Draws in arithmetic progression mod 1 correlate
`1 − 6|d|(1 − |d|)`, so the first draw of replicate `i` correlated 0.95 with that of replicate
`i + 1`. A longer warm-up changes `d` and leaves the progression in place: with 50 warm-up draws
in place of 10, the step is −0.32 and the correlation −0.30. An affine rule for the seeds keeps
the draws in progression too, so no arithmetic scramble helps: a Weyl increment followed by a
Lehmer step left the first draws correlated at −0.27. Matsumoto et al. (2007) describe the same
defect in the initialisation of other linear generators. The first draws of the seeds the loops
use correlate by 0.02 or less, either way. This changed every number pilotr produced before 0.3.

## Sampling algorithms (identical across all implementations)

The rules below turn the random stream into the values a specification describes. Each is written
as the operations an implementation performs, in the order it performs them, since a reordering
lands on a different double or consumes the stream differently. In each rule, `u` is the next
uniform and `z` the next normal, taken from the stream where they appear, and `⌊x⌋` is the
largest whole number not above `x`. Every operation is carried out in double precision and
rounded to nearest.

### Partial crossing

With `n` items, of which each subject sees `m`, step 0 of the draw order samples one subject's
items by a partial Fisher–Yates shuffle, counting positions from 0:

```
pool ← [1, 2, …, n]
for k = 0, 1, …, m − 1:
    j ← k + ⌊u · (n − k)⌋
    swap pool[k] and pool[j]
items ← pool[0 … m − 1], sorted ascending
```

### Blocks of units

A between factor with `L` levels over `N` units gives unit `t`, counted from 1, the level with
index `⌊(t − 1) · L / N⌋`, counted from 0. An additional grouping factor with `K` groups over `N`
units puts unit `t` in group `⌊(t − 1) · K / N⌋`, counted from 0 and written to the data counted
from 1. Both are exact integer arithmetic and draw nothing.

From version 0.4, the `k` between factors of one unit, in listed order with `L_1, …, L_k`
levels, are allocated together. With one factor, the result is the level above.

```
C ← L_1 · L_2 · … · L_k
c ← ⌊(t − 1) · C / N⌋
for j = k, k − 1, …, 1:
    level_j ← c mod L_j
    c ← ⌊c / L_j⌋
```

### Predictor values

A predictor draws its latent value as

```
normal:   value ← mean + sd · z
uniform:  value ← min + (max − min) · u
```

where `mean` defaults to 0 and `sd` to 1. With a `reliability` `ρ` below 1, one further normal per
value gives the observed value. `m` and `s` are the population mean and standard deviation of the
latent variable: `mean` and `sd` for a normal predictor, and `(min + max) / 2` and
`(max − min) / sqrt(12)` for a uniform one.

```
observed ← m + (value − m + s · sqrt((1 − ρ) / ρ) · z) · sqrt(ρ)
```

The latent value enters the linear predictor and the observed one the data (see
[Reliability](#reliability)).

### Random effects

A covariance `Σ` of `q` terms, formed as under [Random effects](#random-effects-per-unit), is
factorised row by row by the Cholesky–Banachiewicz algorithm. Every sum is a left fold taken in
index order, and an empty one is 0.

```
tol ← ε · q · max(Σ[1,1], …, Σ[q,q], 1)          with ε = 2^−52
for i = 1, …, q:
    for j = 1, …, i:
        s ← L[i,1]·L[j,1] + … + L[i,j−1]·L[j,j−1]
        if i = j:
            d ← Σ[i,i] − s
            refuse if d < −tol, otherwise L[i,i] ← sqrt(max(d, 0))
        otherwise:
            L[i,j] ← (Σ[i,j] − s) / L[j,j] if L[j,j] ≠ 0, otherwise 0
```

Each unit, or each group of an additional grouping factor, then takes `b ← L · z` for `q` new
normals, each element of `b` a left fold over all `q` terms of its row.

### Inverse logit

```
invlogit(x) ← 1 / (1 + exp(−x))             if x ≥ 0
invlogit(x) ← e / (1 + e), with e ← exp(x)  if x < 0
```

The exponent is never positive, so `exp()` cannot overflow.

### Responses

With `η` the linear predictor of the row, the response is drawn as below, and a rounded one then
follows [Rounding](#rounding-identical-across-all-implementations).

```
gaussian            y ← η + σ · z
lognormal           y ← exp(η + σ · z)
shifted_lognormal   y ← shift + exp(η + σ · z)
exgaussian          y ← (η + σ · z) − β · (log(u) + 1)        z first, then u
bernoulli           y ← 1 if u < invlogit(η), otherwise 0
poisson             y ← poisson(exp(η), u)
ordinal             y ← the smallest k with u ≤ invlogit(θ_k − η), or K if there is none
beta                y ← beta(η, φ)
```

Ordinal categories are counted from 1, and `K` is one more than the number of thresholds. A
lognormal or shifted lognormal response is refused when its `exp()` overflows. The Poisson and
Beta rules refuse the cases set out under [Response families](#response-families):

```
poisson(λ, u):
    refuse if λ > 708.3964185322641
    p ← exp(−λ);  c ← p;  k ← 0
    while u > c and k < 1,000,000:
        k ← k + 1;  p ← (p · λ) / k;  c ← c + p
    refuse if u > c, otherwise return k

beta(η, φ):
    μ ← invlogit(η);  a ← μ · φ;  b ← (1 − μ) · φ
    refuse if a = 0 or b = 0
    x ← gamma(a);  y ← gamma(b)
    refuse if x + y = 0
    p ← x / (x + y)
    return 1 − 2^−53 if p = 1, 2^−1074 if p = 0, otherwise p
```

The Gamma variates come from the method of Marsaglia and Tsang (2000). Its acceptance test is the
logarithmic one, with no squeeze step, and the right-hand side is evaluated from left to right. A
shape below 1 is boosted, and the boosting uniform is drawn after the Gamma(shape + 1) variate:

```
gamma(a):
    if a < 1:
        g ← gamma(a + 1)
        return g · u^(1/a)
    d ← a − 1/3;  c ← 1 / sqrt(9 · d)
    repeat:
        x ← z;  t ← 1 + c · x;  v ← t · t · t
        if v ≤ 0: begin the next proposal, drawing no uniform
        if log(u) < 0.5 · x · x + d − d · v + d · log(v): return d · v
```

`u^(1/a)` is the maths library's `pow()`, as Python's `**` computes it. R's `^` calls `pow()`
too, with two exceptions. It computes a square as `u · u`, and, from R 4.6.0, it computes a cube
or a fourth power by repeated multiplication, which rounds two or three times (see
[Scope of the guarantee](#scope-of-the-guarantee)).

A sweep runs the same loops at each point of a grid, such as a list of sample sizes, and gives
each point a seed of its own by the same rule. Point `j` of a sweep over `J` values, counting from
1, takes the `j`-th of the seeds `replicate_seeds(seed, J)` returns. Its replicates take their
seeds from that one, so that no two points share their random numbers. The seeds are drawn in
turn, so the first `J` of a longer sweep are the seeds of a shorter one, and appending points to
a grid leaves the earlier points' results unchanged. A point's seed follows its position in the
grid, though, not the value at that position. A sweep over the seed itself sets each point's seed
to the value swept. Sweeps up to R 0.3.1 and Python 0.3.0 gave every point the specification's own
seed, so replicate `i` drew the same numbers at every point, and the points of a curve shared
their noise. The curve solvers fit the points as independent, and with shared seeds their nominal
95% intervals covered the true sample size in 87.8% of solves on a tight grid. Both packages keep
that behaviour as an option, `.seeds = "shared"` in R's `sweep_spec()` and `seeds = "shared"` in
its curve functions and in Python's `power_curve`.

## Scope of the guarantee

Identical data means bit-identical. The guarantee holds wherever R and Python compute `exp()`,
`log()` and `pow()` with one maths library, as on Linux, where continuous integration runs the
parity harness and every case agrees to the last bit. One Beta case, set out below, is the
exception. IEEE-754 requires correct rounding for addition, subtraction, multiplication,
division and square root, so a value built from those alone is the same double everywhere. The
uniform stream, AS 241's central branch and the Cholesky factors are such values, and the
linear predictor is formed from its terms by those operations alone. The standard leaves the
rounding of `exp()`, `log()` and `pow()` to the library, and on Windows R and CPython do not
share one. R takes them from the MinGW-w64 toolchain it is built with, and CPython from
Microsoft's Universal C Runtime. They were compared on Windows 11 with R 4.6.1 and Python 3.13,
on arguments in the ranges pilotr's families use. The two returned different doubles for `exp()`
on 0.47% of 700,000 arguments, for `log()` on 0.13% of 700,000 and for `pow()` on 0.11% of
200,000. None differed by more than one unit in the last place (ulp). macOS has not been
measured.

On Windows, the cause can be isolated. Taking a shifted-lognormal design and switching only its
family to `gaussian`, so that the seed, the random-effect structure, the linear predictor and the
entire draw sequence are unchanged, gives bit-identical output in both languages, while the
lognormal original differs in a handful of rows. There, `exp()` is the only difference.

Those last bits reach the data wherever a value passes through one of the three functions, and
every family has such a path. The figures below were measured on the same machine.

* Normal draws, the Gaussian family's included. AS 241's tail branch calls `log()` for about 15%
  of draws, and 1.4e−5 of normal draws differed, by 1 to 3 ulp. A Gaussian design of 50,000 rows
  differed in one row by 1 ulp, and so did a crossed design of 40,000 rows.
* The lognormal families and the ex-Gaussian, which apply `exp()` or `log()` to every response.
  Of 50,000 rows, 213 differed for the lognormal family, 124 for the shifted lognormal and 3 for
  the ex-Gaussian, each by 1 ulp.
* An unrounded Beta response, which differed by up to 35 ulp in 54 of 40,000 rows with `φ = 0.8`.
  A 1-ulp difference in `μ`, from `exp()` in the inverse logit, reaches the exponent of the Gamma
  sampler's boost `u^(1/shape)`, which magnifies a relative difference in the exponent by the size
  of the result's logarithm.
* The discrete families, which decide each outcome by comparing a uniform with a value built from
  `exp()`. A last-bit difference changes the outcome only when the uniform falls between the two
  versions of that value, and then by a whole category. None of 50,000 rows differed for any of
  the three, but nothing rules it out.

The Beta case can differ on any platform, Linux included. From R 4.6.0, R's `^` computes a cube or
a fourth power by repeated multiplication, rounding two or three times, where Python's `**` calls
`pow()`. On Windows, the two differed for 26% of uniforms raised to the third power and 34%
raised to the fourth. A Beta response is affected wherever a Gamma shape below 1 has a reciprocal
of exactly 3 or 4, so that its boost is a cube or a fourth power. That takes a linear predictor
shared by every row of a cell, as in a design without random effects or continuous predictors.
With `η = 0` and `φ = 0.5`, both shapes are 1/4, and 251 of 1,000 such rows differed on Windows,
by up to 3 ulp. R before 4.6.0 computed those powers with the maths library, so R's own data for
such a design changed with R 4.6.0.

Rounding has removed every continuous difference observed so far, in each rounded case of the
parity harness and in the designs above rounded to four decimals, but it guarantees nothing. The
rounding rule is shared (see [Rounding](#rounding-identical-across-all-implementations)), yet a
last-bit difference that straddles a rounding boundary still rounds two ways. `round` is refused
for the discrete families in any case, where a difference is a whole category.

Agreement to the last bit in every family has one route: run both languages on a platform where
they share a maths library, as on Linux, outside that Beta case. Short of that, the Gaussian
family confines any difference to the last bits of a rare normal draw, where a discrete family
can move a whole category.

The R package's `build_spec()`, which builds the no-code app's specifications, rounds a response
of the Gaussian and the two lognormal families to four decimals and leaves a `beta` response
unrounded. Rounding moves every proportion near 0 or 1 onto the nearest value inside the interval
at that many decimals. At four decimals, 9.8% of the 40,000 rows of the Beta design above landed
on 0.0001 or 0.9999, and rows at such a bound dominate a Beta likelihood (see `?brms_bridge` in
R). Unrounded, the data keep the distribution the specification describes.

Within one language, the guarantee holds on one platform: a given build of R or Python produces
the same data from the same specification and seed every time. Across operating systems, it
holds on the same terms as across languages, since each language takes `exp()`, `log()` and
`pow()` from the platform's maths library. On Windows, CPython's `exp()` returned a double other
than the correctly rounded one for 0.48% of the 700,000 arguments above, and R's for 0.03%. A run
on Windows and a run on Linux can therefore differ in their last bits even within one language.

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
* Marsaglia, G., & Tsang, W. W. (2000). A simple method for generating gamma variables. *ACM
  Transactions on Mathematical Software, 26*(3), 363–372. https://doi.org/10.1145/358407.358414
* Matsumoto, M., Wada, I., Kuramoto, A., & Ashihara, H. (2007). Common defects in initialization
  of pseudorandom number generators. *ACM Transactions on Modeling and Computer Simulation,
  17*(4), Article 15. https://doi.org/10.1145/1276927.1276928
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
