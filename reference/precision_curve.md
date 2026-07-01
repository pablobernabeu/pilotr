# Precision and ROPE curve over sample size

Sweep the number of subjects and report the ROPE decision probabilities
at each size, to identify the minimum analysable *N* at which a focal
effect reaches a determinate decision with a target probability (for
example 0.90). Calls \[precision_design()\] and so requires the \`lme4\`
package.

## Usage

``` r
precision_curve(
  spec,
  focal,
  subject_ns,
  formula = NULL,
  prep = NULL,
  rope = 0.05,
  n_sims = 60
)
```

## Arguments

- spec:

  A design specification (path or list).

- focal:

  A named numeric vector (coefficient name to true value) or character
  vector of focal coefficient names.

- subject_ns:

  A numeric vector of subject counts to evaluate.

- formula:

  Optional \`lme4\` formula; if \`NULL\` it is derived via
  \[model_formula()\].

- prep:

  Optional data-preparation function; if \`NULL\` it is derived via
  \[model_data()\].

- rope:

  Half-width of the region of practical equivalence.

- n_sims:

  Number of Monte Carlo replicates per sample size.

## Value

A data frame with one row per focal effect and sample size, adding an
\`n_subject\` column to the columns returned by \[precision_design()\].

## Examples

``` r
# \donttest{
if (requireNamespace("lme4", quietly = TRUE)) {
  spec <- build_spec(list(name = "pr", seed = 1, design_kind = "within",
    include_items = TRUE, n_subject = 12, n_item = 12, factor_name = "cond",
    lev1 = "a", lev2 = "b", intercept = 6, effect = 0.05,
    subj_int_sd = 0.12, subj_slope_sd = 0.04, subj_corr = 0.2,
    item_int_sd = 0.08, item_slope_sd = 0.02, item_corr = -0.1,
    family = "shifted_lognormal", resp_name = "", sigma = 0.3, shift = 200))
  precision_curve(spec, focal = c(effect = 0.05), subject_ns = c(12, 18), n_sims = 8)
}
#>   n_subject  param true mean_ci_width p_meaningful p_equivalent n_converged
#> 1        12 effect 0.05     0.1642526        0.000            0           8
#> 2        18 effect 0.05     0.1306036        0.125            0           8
# }
```
