---
title: "Getting started with pilotr"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Getting started with pilotr}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




``` r
library(pilotr)
```

## One specification, three interfaces

**pilotr** is built around a single idea: a *portable design specification* that fully
describes a data-generating process and is consumed identically by three interfaces — a
no-code web app, this R package, and a Python package. The specification is the unit of
reproducibility: not a screenshot of a GUI, nor a language-specific script, but a complete,
executable description you can version-control, attach to a preregistration, share with a
collaborator, and run unchanged in R or Python to get **bit-identical** data.

This vignette walks through the core loop: **describe → simulate → inspect → export**. The
other vignettes go deeper into [power](power-analysis.html),
[precision/ROPE design analysis](precision-rope.html),
[cross-language reproducibility](cross-language.html), and the
[response families](response-families.html).

## Describe a design

The specification is a small list (serialised to JSON). You can write it by hand, point and
click in the web app, or build it from a flat set of inputs with `build_spec()`. Here is a
crossed psycholinguistic design: subjects respond to items in a *related* or *unrelated*
condition, with reaction time as the outcome, crossed by-subject and by-item random
intercepts and slopes, and a small priming effect on the log scale.


``` r
spec <- build_spec(list(
  name = "priming", seed = 2024,
  design_kind = "within", include_items = TRUE,
  n_subject = 24, n_item = 20,
  factor_name = "condition", lev1 = "related", lev2 = "unrelated",
  intercept = 6.0, effect = 0.05,
  subj_int_sd = 0.12, subj_slope_sd = 0.04, subj_corr = 0.2,
  item_int_sd = 0.08, item_slope_sd = 0.02, item_corr = -0.1,
  family = "shifted_lognormal", resp_name = "RT", sigma = 0.30, shift = 200
))
```

The portable artifact is plain JSON — this is exactly what the web app downloads and what
the Python package reads:


``` r
cat(spec_json(spec))
#> {
#>   "name": "priming",
#>   "seed": 2024,
#>   "units": {
#>     "subject": {
#>       "n": 24
#>     },
#>     "item": {
#>       "n": 20
#>     }
#>   },
#>   "factors": [
#>     {
#>       "name": "condition",
#>       "levels": ["related", "unrelated"],
#>       "contrasts": {
#>         "effect": [-0.5, 0.5]
#>       },
#>       "vary_within": ["subject", "item"]
#>     }
#>   ],
#>   "fixed": {
#>     "intercept": 6,
#>     "coefficients": {
#>       "effect": 0.05
#>     }
#>   },
#>   "random": {
#>     "subject": {
#>       "intercept_sd": 0.12,
#>       "slopes": {
#>         "effect": 0.04
#>       },
#>       "correlations": {
#>         "intercept,effect": 0.2
#>       }
#>     },
#>     "item": {
#>       "intercept_sd": 0.08,
#>       "slopes": {
#>         "effect": 0.02
#>       },
#>       "correlations": {
#>         "intercept,effect": -0.1
#>       }
#>     }
#>   },
#>   "response": {
#>     "family": "shifted_lognormal",
#>     "name": "RT",
#>     "sigma": 0.3,
#>     "round": 4,
#>     "shift": 200
#>   }
#> }
```

You can read a specification back from a file (or string) with `load_spec()`, so a design
authored anywhere runs everywhere.

## Simulate

`simulate_design()` turns a specification (a list **or** a path to a JSON file) into an
analysis-ready data frame:


``` r
data <- simulate_design(spec)
nrow(data)
#> [1] 960
head(data)
#>   subject item condition       RT
#> 1       1    1   related 431.0559
#> 2       1    1 unrelated 488.0423
#> 3       1    2   related 665.4629
#> 4       1    2 unrelated 771.0603
#> 5       1    3   related 480.7395
#> 6       1    3 unrelated 877.9173
```

The design is crossed, so every subject sees every item in both conditions
(960 rows here).

## Reproducible by construction

pilotr ships a shared random-number generator, so the *same specification and seed* produce
the *same data every time* — in this R session, in a fresh one, in Python, and on any
machine:


``` r
isTRUE(all.equal(simulate_design(spec), simulate_design(spec)))
#> [1] TRUE
```

Change the seed and the data change; change nothing and they never do. See the
[cross-language vignette](cross-language.html) for how this holds across R and Python.

## Export a self-contained script

`generate_r_script()` emits a stand-alone R script that embeds the specification and
reproduces the design with no external files — handy for an appendix or a preregistration:


``` r
cat(generate_r_script(spec))
#> # Reproducible simulation exported by pilotr.
#> # install.packages("pilotr")   # once available; then run this script as-is.
#> library(pilotr)
#> 
#> spec <- list(name = "priming", seed = 2024L, units = list(subject = list(
#>     n = 24L), item = list(n = 20L)), factors = list(list(name = "condition", 
#>     levels = c("related", "unrelated"), contrasts = list(effect = c(-0.5, 
#>     0.5)), vary_within = c("subject", "item"))), fixed = list(
#>     intercept = 6, coefficients = list(effect = 0.05)), random = list(
#>     subject = list(intercept_sd = 0.12, slopes = list(effect = 0.04), 
#>         correlations = list("intercept,effect" = 0.2)), item = list(
#>         intercept_sd = 0.08, slopes = list(effect = 0.02), correlations = list(
#>             "intercept,effect" = -0.1))), response = list(family = "shifted_lognormal", 
#>     name = "RT", sigma = 0.3, round = 4L, shift = 200))
#> 
#> data <- simulate_design(spec)              # analysis-ready data frame
#> # write.csv(data, "data.csv", row.names = FALSE)
#> # pow  <- power_mixed(spec, n_sims = 200)   # simulation-based power + Type S/M
```

## Where to go next

- **[Power & design analysis](power-analysis.html)** — simulation-based power with Type S and
  Type M errors, for two-group and crossed mixed-effects designs.
- **[Precision / ROPE design analysis](precision-rope.html)** — plan by precision against a
  region of practical equivalence, and sweep sample size.
- **[Cross-language reproducibility](cross-language.html)** — the shared generator that makes
  R and Python agree to the last bit.
- **[Response families](response-families.html)** — Gaussian, reaction-time, accuracy, count,
  ordinal, and proportion outcomes.
