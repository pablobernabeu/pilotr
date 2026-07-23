# pilotr (R)

A twin of the [Python package of the same
name](https://pablobernabeu.github.io/pilotr/py/): the two share the
design specification and the random-number generator, so the same
specification and seed produce bit-identical data in either language.
The R package additionally offers precision/ROPE design analysis, the
`lme4` reference backend for mixed-effects power, the `brms` bridge and
the app launcher.

Simulate experimental and behavioural data from a portable design
specification, with integrated simulation-based power and design
analysis (the Type S and Type M errors, and precision against a region
of practical equivalence).

pilotr lets you pilot a study before you run it. Describe the design you
plan to collect, with its groups, conditions, sample sizes, effect sizes
and outcome family, and pilotr generates the data that design would
produce. You can then check the study’s power and design analysis before
gathering anything. A single specification drives three interchangeable
interfaces, namely a no-code web app, this R package and the Python
package.

## Installation

``` r

# install.packages("remotes")
remotes::install_github("pablobernabeu/pilotr", subdir = "r/pilotr")
```

## Quick start

The package ships a specification per design family, so the example
below runs as it stands.
[`pilotr_example()`](https://pablobernabeu.github.io/pilotr/reference/pilotr_example.md)
returns the path to one of them, here a crossed by-subject and by-item
reaction-time design.

``` r

library(pilotr)

spec <- load_spec(pilotr_example("crossed_mixed_rt"))
data <- simulate_design(spec)
head(data)
#>   subject item condition       RT
#> 1       1    1   related 668.6564
#> 2       1    1 unrelated 727.3870
#> 3       1    2   related 699.1907
#> 4       1    2 unrelated 502.9760
#> 5       1    3   related 661.0859
#> 6       1    3 unrelated 516.6668
```

The specification carries a seed, so those rows are the same on every
machine and in the Python twin. Point
[`load_spec()`](https://pablobernabeu.github.io/pilotr/reference/load_spec.md)
at your own `design.json` to simulate a design of your own.

Power and design analysis run from the same object. The call below is
not evaluated here because it needs `lme4`, a suggested rather than a
required dependency, and a few hundred model fits.

``` r

power_mixed(spec, n_sims = 200)   # crossed-LMM power + Type S/M (lme4)
cat(generate_r_script(spec))      # a self-contained, reproducible script
```

The [power
article](https://pablobernabeu.github.io/pilotr/articles/power-analysis.html)
shows its output, alongside a power curve over sample size.

## Try it without installing

A serverless build runs entirely in your browser, with no installation
required and no data uploaded. It is available as a [no-code
app](https://pablobernabeu.github.io/pilotr/app/).

## Learn more

The [Get started
article](https://pablobernabeu.github.io/pilotr/articles/getting-started.html)
walks through the core loop of describe, simulate, inspect and export,
and the other articles each go deeper into one part of the workflow. The
full repository, including the specification format and the no-code app,
is at <https://github.com/pablobernabeu/pilotr>.
