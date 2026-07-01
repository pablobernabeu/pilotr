# Build a design specification from a flat list of design inputs

Assemble a portable design specification (a plain list, serialisable
with \[spec_json()\]) from the flat set of inputs collected by the
no-code application: sample sizes, the two-level factor and its levels,
the fixed intercept and effect, the random-effect standard deviations
for within-subject and crossed designs, and the response family with its
parameters.

## Usage

``` r
build_spec(p)
```

## Arguments

- p:

  A named list of design inputs. Common fields are \`name\`, \`seed\`,
  \`n_subject\`, \`design_kind\` (\`"between"\` or \`"within"\`),
  \`include_items\`, \`n_item\`, \`factor_name\`, \`lev1\`, \`lev2\`,
  \`intercept\`, \`effect\`, \`family\`, \`resp_name\`, and family
  parameters such as \`sigma\`, \`shift\`, \`thresholds\`, or \`phi\`;
  within-design random effects use \`subj_int_sd\`, \`subj_slope_sd\`,
  \`subj_corr\`, \`item_int_sd\`, \`item_slope_sd\`, and \`item_corr\`.

## Value

A design specification as a nested list, ready for
\[simulate_design()\], \[spec_json()\], or the power and precision
functions.

## Examples

``` r
build_spec(list(name = "demo", seed = 1, design_kind = "between",
  factor_name = "group", lev1 = "control", lev2 = "treatment", n_subject = 40,
  intercept = 100, effect = 5, family = "gaussian", resp_name = "", sigma = 10))
#> $name
#> [1] "demo"
#> 
#> $seed
#> [1] 1
#> 
#> $units
#> $units$subject
#> $units$subject$n
#> [1] 40
#> 
#> 
#> 
#> $factors
#> $factors[[1]]
#> $factors[[1]]$name
#> [1] "group"
#> 
#> $factors[[1]]$levels
#> [1] "control"   "treatment"
#> 
#> $factors[[1]]$contrasts
#> $factors[[1]]$contrasts$effect
#> [1] -0.5  0.5
#> 
#> 
#> $factors[[1]]$between
#> [1] "subject"
#> 
#> 
#> 
#> $fixed
#> $fixed$intercept
#> [1] 100
#> 
#> $fixed$coefficients
#> $fixed$coefficients$effect
#> [1] 5
#> 
#> 
#> 
#> $response
#> $response$family
#> [1] "gaussian"
#> 
#> $response$name
#> [1] "score"
#> 
#> $response$sigma
#> [1] 10
#> 
#> $response$round
#> [1] 4
#> 
#> 
```
