# Getting started with pilotr

``` r

library(pilotr)
```

## One specification, three interfaces

pilotr is built around a portable design specification that fully
describes a data-generating process and is consumed identically by three
interfaces, namely a no-code web application, this R package and a
Python package. The specification is the unit of reproducibility. It is
a complete, executable description rather than a screenshot of a
graphical interface or a language-specific script. You can place it
under version control, attach it to a preregistration, share it with a
collaborator, and run it unchanged in R or Python to obtain
bit-identical data.

This vignette walks through the core loop of describe, simulate, inspect
and export. The other vignettes, listed at the end, each go deeper into
one part of the workflow.

## Describe a design

The specification is a small list serialised to JSON. You can write it
by hand, point and click in the web application, or build it from a flat
set of inputs with
[`build_spec()`](https://pablobernabeu.github.io/pilotr/reference/build_spec.md).
The example below is a crossed psycholinguistic design. Subjects respond
to items in a related or unrelated condition, with reaction time as the
outcome, crossed by-subject and by-item random intercepts and slopes,
and a small priming effect on the log scale.

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

The portable artefact is plain JSON. This is exactly what the web
application downloads and what the Python package reads.

``` r
cat(spec_json(spec))
{
  "name": "priming",
  "seed": 2024,
  "units": {
    "subject": {
      "n": 24
    },
    "item": {
      "n": 20
    }
  },
  "factors": [
    {
      "name": "condition",
      "levels": ["related", "unrelated"],
      "contrasts": {
        "effect": [-0.5, 0.5]
      },
      "vary_within": ["subject", "item"]
    }
  ],
  "fixed": {
    "intercept": 6,
    "coefficients": {
      "effect": 0.05
    }
  },
  "random": {
    "subject": {
      "intercept_sd": 0.12,
      "slopes": {
        "effect": 0.04
      },
      "correlations": {
        "intercept,effect": 0.2
      }
    },
    "item": {
      "intercept_sd": 0.08,
      "slopes": {
        "effect": 0.02
      },
      "correlations": {
        "intercept,effect": -0.1
      }
    }
  },
  "response": {
    "family": "shifted_lognormal",
    "name": "RT",
    "sigma": 0.3,
    "round": 4,
    "shift": 200
  }
}
```

You can read a specification back from a file (or string) with
[`load_spec()`](https://pablobernabeu.github.io/pilotr/reference/load_spec.md),
so a design authored anywhere runs everywhere. The full specification
format, including every field and the cross-language RNG draw order, is
documented in
[`spec/SPEC.md`](https://github.com/pablobernabeu/pilotr/blob/main/spec/SPEC.md).

## Simulate

[`simulate_design()`](https://pablobernabeu.github.io/pilotr/reference/simulate_design.md)
turns a specification, supplied either as a list or as a path to a JSON
file, into an analysis-ready data frame.

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

The design is crossed, so every subject sees every item in both
conditions (960 rows here).

## Reproducible by construction

pilotr ships a shared random-number generator, so the same specification
and seed produce the same data every time. This holds in the current R
session, in a fresh one, in Python, and on any machine.

``` r

isTRUE(all.equal(simulate_design(spec), simulate_design(spec)))
#> [1] TRUE
```

Changing the seed changes the data, and leaving the specification
unchanged leaves the data unchanged. The [cross-language
vignette](https://pablobernabeu.github.io/pilotr/articles/cross-language.md)
describes how this holds across R and Python.

## Export a self-contained script

[`generate_r_script()`](https://pablobernabeu.github.io/pilotr/reference/generate_r_script.md)
emits a stand-alone R script that embeds the specification and
reproduces the design with no external files. This is convenient for an
appendix or a preregistration.

``` r

cat(generate_r_script(spec))
# Reproducible simulation exported by pilotr.
# install.packages("pilotr")   # once available; then run this script as-is.
library(pilotr)

spec <- list(name = "priming", seed = 2024L, units = list(subject = list(
    n = 24L), item = list(n = 20L)), factors = list(list(name = "condition", 
    levels = c("related", "unrelated"), contrasts = list(effect = c(-0.5, 
    0.5)), vary_within = c("subject", "item"))), fixed = list(
    intercept = 6, coefficients = list(effect = 0.05)), random = list(
    subject = list(intercept_sd = 0.12, slopes = list(effect = 0.04), 
        correlations = list("intercept,effect" = 0.2)), item = list(
        intercept_sd = 0.08, slopes = list(effect = 0.02), correlations = list(
            "intercept,effect" = -0.1))), response = list(family = "shifted_lognormal", 
    name = "RT", sigma = 0.3, round = 4L, shift = 200))

data <- simulate_design(spec)              # analysis-ready data frame
# write.csv(data, "data.csv", row.names = FALSE)
# pow  <- power_mixed(spec, n_sims = 200)   # simulation-based power + Type S/M
```

## Where to go next

- [Response
  families](https://pablobernabeu.github.io/pilotr/articles/response-families.md),
  covering Gaussian, reaction-time, accuracy, count, ordinal and
  proportion outcomes.
- [Power and design
  analysis](https://pablobernabeu.github.io/pilotr/articles/power-analysis.md),
  covering simulation-based power with Type S and Type M errors, for
  two-group and crossed mixed-effects designs.
- [Precision / ROPE design
  analysis](https://pablobernabeu.github.io/pilotr/articles/precision-rope.md),
  which plans by precision against a region of practical equivalence and
  sweeps sample size.
- [Cross-language
  reproducibility](https://pablobernabeu.github.io/pilotr/articles/cross-language.md),
  describing the shared generator that makes R and Python agree to the
  last bit.
- [The no-code
  app](https://pablobernabeu.github.io/pilotr/articles/the-no-code-app.md),
  the point-and-click interface over the same specification, running in
  the browser or locally with
  [`run_app()`](https://pablobernabeu.github.io/pilotr/reference/run_app.md).
