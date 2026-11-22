# Getting started with pilotr

``` r

library(pilotr)
```

## One specification, three interfaces

pilotr is built around a portable design specification that fully
describes a data-generating process. Three interfaces read it
identically, namely a no-code web application, this R package and a
Python package. The specification, an executable description of the
design, is the unit of reproducibility. You can place it under version
control, attach it to a preregistration, share it with a collaborator
and run it unchanged in R or Python to obtain identical data. The data
agree to the last bit wherever R and Python share a maths library, as on
Linux, except in a rare Beta case under R 4.6.0 or later that the
[cross-language
vignette](https://pablobernabeu.github.io/pilotr/r/articles/cross-language.md)
describes. Elsewhere, a value that passes through
[`exp()`](https://rdrr.io/r/base/Log.html),
[`log()`](https://rdrr.io/r/base/Log.html) or `pow()` can differ in its
last bits, which the same vignette quantifies.

This vignette covers the core workflow of describing, simulating,
inspecting and exporting a design. Each of the other vignettes, listed
at the end, covers one part of the workflow in more detail.

## Describe a design

The specification is a small list serialised to JSON. You can write it
by hand, point and click in the web application or build it from a flat
set of inputs with
[`build_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/build_spec.md).
The example below is a crossed psycholinguistic design in which subjects
respond to items in a related or unrelated condition, with reaction time
as the outcome. It has a small priming effect on the log scale and
crossed by-subject and by-item random intercepts and slopes.

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

The JSON form of the specification is what the web application downloads
and what the Python package reads.

``` r

cat(spec_json(spec))
```

    {
      "spec_version": "0.3",
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

You can read a specification back from a file or a string with
[`load_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/load_spec.md).
The specification format, including every field and the cross-language
RNG draw order, is documented in
[`spec/SPEC.md`](https://github.com/pablobernabeu/pilotr/blob/main/spec/SPEC.md).

## Simulate

[`simulate_design()`](https://pablobernabeu.github.io/pilotr/r/reference/simulate_design.md)
turns a specification, supplied either as a list or as a path to a JSON
file, into an analysis-ready data frame.

``` r

data <- simulate_design(spec)
nrow(data)
```

    [1] 960

``` r

head(data)
```

      subject item condition       RT
    1       1    1   related 431.0559
    2       1    1 unrelated 488.0423
    3       1    2   related 665.4629
    4       1    2 unrelated 771.0603
    5       1    3   related 480.7395
    6       1    3 unrelated 877.9173

The design is crossed, so every subject sees every item in both
conditions (960 rows here).

## Reproducible by construction

pilotr provides a shared random-number generator, so the same
specification and seed produce the same data every time. This holds in
the current R session and in a fresh one on the same machine. In Python,
and on another machine, it holds to the last bit wherever the maths
library is the same, except in the rare Beta case under R 4.6.0 or
later.

``` r

identical(simulate_design(spec), simulate_design(spec))
```

    [1] TRUE

Changing the seed changes the data. The [cross-language
vignette](https://pablobernabeu.github.io/pilotr/r/articles/cross-language.md)
describes how reproducibility holds across R and Python.

## Export a self-contained script

[`generate_r_script()`](https://pablobernabeu.github.io/pilotr/r/reference/generate_r_script.md)
emits a stand-alone R script suited to an appendix or a preregistration.
The script reproduces the design with no external files and embeds the
specification as the JSON shown above, which
[`load_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/load_spec.md)
reads back to exactly the same numbers on any build of R.

``` r

cat(generate_r_script(spec))
```

    # Reproducible simulation exported by pilotr.
    # install.packages(c("pilotr", "lme4", "lmerTest"))   # from CRAN; lme4 and lmerTest only for power_mixed()
    library(pilotr)

    # The specification as spec_json() writes it. load_spec() reads it with jsonlite, which
    # recovers every number exactly on every build of R.
    spec <- load_spec(r"(
    {
      "spec_version": "0.3",
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
    )")

    data <- simulate_design(spec)              # analysis-ready data frame
    # write.csv(data, "data.csv", row.names = FALSE)
    # pow <- power_mixed(spec, n_sims = 200)   # needs lme4 and lmerTest

## Where to go next

- [Response
  families](https://pablobernabeu.github.io/pilotr/r/articles/response-families.md),
  covering Gaussian, reaction-time, accuracy, count, ordinal and
  proportion outcomes.
- [Worked
  examples](https://pablobernabeu.github.io/pilotr/r/articles/worked-examples.md),
  which simulates each specification bundled with the package.
- [Power and design
  analysis](https://pablobernabeu.github.io/pilotr/r/articles/power-analysis.md),
  covering simulation-based power with Type S and Type M errors, for
  two-group and crossed mixed-effects designs.
- [Precision / ROPE design
  analysis](https://pablobernabeu.github.io/pilotr/r/articles/precision-rope.md),
  which plans by precision against a region of practical equivalence and
  sweeps sample size.
- [Cross-language
  reproducibility](https://pablobernabeu.github.io/pilotr/r/articles/cross-language.md),
  describing the shared generator that makes R and Python agree, and
  where their data can differ in the last bits.
- [The no-code
  app](https://pablobernabeu.github.io/pilotr/r/articles/the-no-code-app.md),
  the point-and-click interface over the same specification, running in
  the browser or locally with
  [`run_app()`](https://pablobernabeu.github.io/pilotr/r/reference/run_app.md).
