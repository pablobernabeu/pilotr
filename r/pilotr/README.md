# pilotr <img src="man/figures/logo.png" align="right" height="139" alt="pilotr hex logo" />

Simulate experimental and behavioural data from a portable design specification, with
integrated simulation-based power and design analysis (the Type S and Type M errors, and
precision against a region of practical equivalence).

pilotr lets you pilot a study before you run it. Describe the design you plan to collect, with
its groups, conditions, sample sizes, effect sizes and outcome family, and pilotr generates the
data that design would produce. You can then check the study's power and design analysis before
gathering anything. A single specification drives three interchangeable interfaces, namely a
no-code web app, this R package and a Python package; given the same specification and seed, all
three produce bit-identical data.

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
uploaded. It is available as a [live demo](https://pablobernabeu.github.io/pilotr/demo/).

## Learn more

This is the R package's documentation. Its companion
**[Python package](https://pablobernabeu.github.io/pilotr/py/)** mirrors it exactly and produces
identical data from the same specification. The full repository, including the specification
format and the no-code app, is at <https://github.com/pablobernabeu/pilotr>.
