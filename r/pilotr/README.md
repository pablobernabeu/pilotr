# pilotr <small>(R)</small> <img src="man/figures/logo.png" align="right" height="139" alt="pilotr hex logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/pablobernabeu/pilotr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pablobernabeu/pilotr/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/license/MIT)
<!-- badges: end -->

A twin of the [Python package of the same name](https://pablobernabeu.github.io/pilotr/py/): the two share the design specification and the random-number generator, so the same specification and seed produce bit-identical data in either language. The R package additionally offers precision/ROPE design analysis, the `lme4` reference backend for mixed-effects power, the `brms` bridge and the app launcher.

Simulate experimental and behavioural data from a portable design specification, with
integrated simulation-based power and design analysis (the Type S and Type M errors, and
precision against a region of practical equivalence).

pilotr lets you pilot a study before you run it. Describe the design you plan to collect, with
its groups, conditions, sample sizes, effect sizes and outcome family, and pilotr generates the
data that design would produce. You can then check the study's power and design analysis before
gathering anything. A single specification drives three interchangeable interfaces, namely a
no-code web app, this R package and the Python package.

## Installation

```r
# install.packages("remotes")
remotes::install_github("pablobernabeu/pilotr", subdir = "r/pilotr")
```

## Quick start

```r
library(pilotr)

spec <- load_spec("design.json")          # a pilotr JSON design spec
data <- simulate_design(spec)             # analysis-ready data frame
pow  <- power_mixed(spec, n_sims = 200)   # crossed-LMM power + Type S/M (lme4)
cat(generate_r_script(spec))              # a self-contained, reproducible script
```

## Try it without installing

A serverless build runs entirely in your browser, with no installation required and no data
uploaded. It is available as a [no-code app](https://pablobernabeu.github.io/pilotr/app/).

## Learn more

The [Get started article](https://pablobernabeu.github.io/pilotr/articles/getting-started.html)
walks through the core loop of describe, simulate, inspect and export, and the other articles
each go deeper into one part of the workflow. The full repository, including the specification
format and the no-code app, is at <https://github.com/pablobernabeu/pilotr>.
