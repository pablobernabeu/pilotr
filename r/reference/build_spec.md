# Build a design specification from a flat list of design inputs

Assemble a portable design specification (a plain list, serialisable
with
[`spec_json()`](https://pablobernabeu.github.io/pilotr/r/reference/spec_json.md))
from the flat set of inputs collected by the no-code application: sample
sizes, the two-level factor and its levels, the fixed intercept and
effect, the random-effect standard deviations for within-subject and
crossed designs, and the response family with its parameters.

## Usage

``` r
build_spec(p)
```

## Arguments

- p:

  A named list of design inputs. Common fields are `name`, `seed`,
  `n_subject`, `design_kind` (`"between"` or `"within"`),
  `include_items`, `n_item`, `factor_name`, `lev1`, `lev2`, `intercept`,
  `effect`, `family`, `resp_name`, and family parameters such as
  `sigma`, `shift`, `beta`, `thresholds`, or `phi`; within-design random
  effects use `subj_int_sd`, `subj_slope_sd`, `subj_corr`,
  `item_int_sd`, `item_slope_sd`, and `item_corr`.
  `counterbalance = TRUE` counterbalances the factor of a within design
  with items, and is ignored in any other design.

## Value

A design specification as a nested list, ready for
[`simulate_design()`](https://pablobernabeu.github.io/pilotr/r/reference/simulate_design.md),
[`spec_json()`](https://pablobernabeu.github.io/pilotr/r/reference/spec_json.md),
or the power and precision functions.

## Details

A response of the Gaussian, lognormal, shifted lognormal or ex-Gaussian
family is rounded to four decimals, and a `beta` response is left
unrounded. Rounding moves every proportion near 0 or 1 onto the nearest
value inside the interval at that many decimals. At four decimals, 9.8%
of the 40,000 rows of a design with `phi = 0.8` landed on 0.0001 or
0.9999, and rows at such a bound dominate a Beta likelihood (see
[`brms_bridge()`](https://pablobernabeu.github.io/pilotr/r/reference/brms_bridge.md)).
Unrounded, the data keep the distribution the specification describes.
They can then differ between R and 'Python' in their last bits where the
two do not share a maths library, as the cross-language vignette sets
out.

Counterbalanced lists show each subject each item at one level of the
factor, rotating the levels over subjects and items as a Latin square
does, so a design has half the rows of the fully crossed one. They set
`counterbalance` on the factor, a field of specification version 0.4
described under Counterbalanced lists in the specification.

The specification declares the lowest `spec_version` that reads it. Only
counterbalanced lists need a version above 0.3, so the releases of
either package that read 0.3 read every other design built here.

## Examples

``` r
build_spec(list(name = "demo", seed = 1, design_kind = "between",
  factor_name = "group", lev1 = "control", lev2 = "treatment", n_subject = 40,
  intercept = 100, effect = 5, family = "gaussian", resp_name = "", sigma = 10))
#> $spec_version
#> [1] "0.3"
#> 
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
