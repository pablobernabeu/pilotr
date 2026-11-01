# Derive the lmer formula implied by a specification

Construct the mixed-model formula (with response `.y`, as produced by
[`model_data()`](https://pablobernabeu.github.io/pilotr/r/reference/model_data.md))
from the fixed-effect coefficients and random-effects structure of a
specification. Interaction coefficients written `a:b` become formula
terms `a_b`.

## Usage

``` r
model_formula(spec, structure = c("maximal", "zcp", "intercepts"))
```

## Arguments

- spec:

  A design specification (path or list). A path is read with
  [`load_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/load_spec.md),
  which validates it.

- structure:

  The random-effects structure: `"maximal"`, the default, `"zcp"` or
  `"intercepts"` (see Details).

## Value

A [stats::formula](https://rdrr.io/r/stats/formula.html) object suitable
for fitting with `lme4` or `lmerTest`.

## Details

The three structures are the usual steps of a rule that simplifies a
model when its fit is singular or fails to converge (Barr et al., 2013;
Bates et al., 2015; Matuschek et al., 2017). `"maximal"` has every
random slope the specification declares, correlated with the intercept
where the specification correlates them, and uncorrelated (`||`) where
it does not. `"zcp"`, the zero-correlation-parameter model, keeps those
slopes and estimates no correlation among them, writing `||` for every
grouping factor with a slope. Since
[`model_data()`](https://pablobernabeu.github.io/pilotr/r/reference/model_data.md)
supplies numeric contrast columns, `||` does remove the correlations,
which it does not for a factor. A grouping factor without slopes keeps
`(1 | g)`, since lme4 refuses `(1 || g)`. `"intercepts"` keeps a random
intercept for each grouping factor and drops every slope.
[`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md)
and
[`precision_design()`](https://pablobernabeu.github.io/pilotr/r/reference/precision_design.md)
take these names as the rungs of such a rule.

A specification with neither a fixed coefficient nor a random effect
implies no model, and is refused in the words the 'Python' twin uses.

## References

Barr, D. J., Levy, R., Scheepers, C. and Tily, H. J. (2013). Random
effects structure for confirmatory hypothesis testing: Keep it maximal.
*Journal of Memory and Language*, 68(3), 255-278.
[doi:10.1016/j.jml.2012.11.001](https://doi.org/10.1016/j.jml.2012.11.001)

Bates, D., Kliegl, R., Vasishth, S. and Baayen, H. (2015). Parsimonious
mixed models. *arXiv*.
[doi:10.48550/arXiv.1506.04967](https://doi.org/10.48550/arXiv.1506.04967)

Matuschek, H., Kliegl, R., Vasishth, S., Baayen, H. and Bates, D.
(2017). Balancing Type I error and power in linear mixed models.
*Journal of Memory and Language*, 94, 305-315.
[doi:10.1016/j.jml.2017.01.001](https://doi.org/10.1016/j.jml.2017.01.001)

## Examples

``` r
spec <- build_spec(list(name = "d", seed = 1, design_kind = "within",
  include_items = TRUE, n_subject = 10, n_item = 8, factor_name = "cond",
  lev1 = "a", lev2 = "b", intercept = 6, effect = 0.05,
  subj_int_sd = 0.12, subj_slope_sd = 0.04, subj_corr = 0.2,
  item_int_sd = 0.08, item_slope_sd = 0.02, item_corr = -0.1,
  family = "shifted_lognormal", resp_name = "", sigma = 0.3, shift = 200))
model_formula(spec)
#> .y ~ effect + (1 + effect | subject) + (1 + effect | item)
model_formula(spec, "zcp")
#> .y ~ effect + (1 + effect || subject) + (1 + effect || item)
model_formula(spec, "intercepts")
#> .y ~ effect + (1 | subject) + (1 | item)
```
