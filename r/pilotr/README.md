# pilotr

Simulate experimental and behavioral data from a portable design specification, with
integrated simulation-based power and design analysis (Type S / Type M errors, and
precision against a region of practical equivalence). A single specification drives three
interchangeable interfaces, namely a no-code web app, this R package, and a Python
package. Given the same specification and seed, all three produce bit-identical data.

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

A serverless build runs entirely in your browser, with no installation required and no
data uploaded. It is available as a [live demo](https://pablobernabeu.github.io/pilotr/demo/).

## Learn more

The full toolkit, comprising the specification format, the Python package, the no-code app,
the validation suite, and the manuscript, is available at
<https://github.com/pablobernabeu/pilotr>.
