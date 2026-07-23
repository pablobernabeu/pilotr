# Worked examples

``` r

library(pilotr)
```

pilotr ships one ready-to-run specification per design family. The same
JSON files drive the Python twin and the no-code app, so a design
authored once runs unchanged across all three.
[`pilotr_example()`](https://pablobernabeu.github.io/pilotr/reference/pilotr_example.md)
lists them, and returns the path to each for
[`load_spec()`](https://pablobernabeu.github.io/pilotr/reference/load_spec.md).

``` r

pilotr_example()
#> [1] "beta_proportion"         "between_2group_gaussian"
#> [3] "crossed_mixed_rt"        "nested_clusters"        
#> [5] "ordinal_likert_between"  "partial_crossing"       
#> [7] "poisson_counts_between"  "reading_time_continuous"
```

Each specification is simulated below, showing the family it draws from
and the first rows of the data it produces. The specification format
itself is covered in the [Get
started](https://pablobernabeu.github.io/pilotr/articles/getting-started.md)
article.

``` r

desc <- c(
  between_2group_gaussian = "Two-group between-subjects Gaussian.",
  crossed_mixed_rt        = "Crossed by-subject and by-item reaction times, shifted lognormal.",
  beta_proportion         = "Bounded proportions through the Beta family.",
  ordinal_likert_between  = "Five-point Likert responses via a cumulative-logit model.",
  poisson_counts_between  = "Count outcomes through a log link.",
  reading_time_continuous = "A continuous predictor with a lognormal reading-time outcome.",
  nested_clusters         = "Subjects nested in higher-level clusters, an extra grouping factor.",
  partial_crossing        = "Each subject sees a sampled subset of items."
)

for (name in pilotr_example()) {
  spec <- load_spec(pilotr_example(name))
  d <- simulate_design(spec)
  blurb <- if (name %in% names(desc)) desc[[name]] else ""
  cat(sprintf("\n### %s\n\n", name))
  cat(sprintf("%s The `%s` family, %d rows.\n\n",
              blurb, spec$response$family, nrow(d)))
  # These tables are printed data, so they take the site's code size rather than the prose
  # size a pipe table would inherit. The class is what the stylesheet keys on, and
  # table.attr reaches the output only for format = "html"; pipe output discards it.
  cat(knitr::kable(head(d, 4), format = "html",
                   table.attr = 'class="table data-output"'), sep = "\n")
  cat("\n\n")
}
```

### beta_proportion

Bounded proportions through the Beta family. The `beta` family, 200
rows.

| subject | group   |    prop |
|--------:|:--------|--------:|
|       1 | control | 0.29928 |
|       2 | control | 0.54044 |
|       3 | control | 0.25637 |
|       4 | control | 0.33684 |

### between_2group_gaussian

Two-group between-subjects Gaussian. The `gaussian` family, 64 rows.

| subject | group   |    score |
|--------:|:--------|---------:|
|       1 | control |  95.7157 |
|       2 | control |  90.0921 |
|       3 | control | 119.1958 |
|       4 | control |  86.4388 |

### crossed_mixed_rt

Crossed by-subject and by-item reaction times, shifted lognormal. The
`shifted_lognormal` family, 1440 rows.

| subject | item | condition |       RT |
|--------:|-----:|:----------|---------:|
|       1 |    1 | related   | 668.6564 |
|       1 |    1 | unrelated | 727.3870 |
|       1 |    2 | related   | 699.1907 |
|       1 |    2 | unrelated | 502.9760 |

### nested_clusters

Subjects nested in higher-level clusters, an extra grouping factor. The
`gaussian` family, 4800 rows.

| subject | item | site | condition |       y |
|--------:|-----:|-----:|:----------|--------:|
|       1 |    1 |    1 | a         | 0.16487 |
|       1 |    1 |    1 | b         | 2.11132 |
|       1 |    2 |    1 | a         | 0.72468 |
|       1 |    2 |    1 | b         | 1.70270 |

### ordinal_likert_between

Five-point Likert responses via a cumulative-logit model. The `ordinal`
family, 400 rows.

| subject | group   | rating |
|--------:|:--------|-------:|
|       1 | control |      2 |
|       2 | control |      4 |
|       3 | control |      3 |
|       4 | control |      3 |

### partial_crossing

Each subject sees a sampled subset of items. The `gaussian` family, 1440
rows.

| subject | item | condition |       rt |
|--------:|-----:|:----------|---------:|
|       1 |    3 | a         | 393.4013 |
|       1 |    3 | b         | 532.2501 |
|       1 |    5 | a         | 524.1783 |
|       1 |    5 | b         | 632.9765 |

### poisson_counts_between

Count outcomes through a log link. The `poisson` family, 2000 rows.

| subject | group   | count |
|--------:|:--------|------:|
|       1 | control |     3 |
|       2 | control |     5 |
|       3 | control |     3 |
|       4 | control |     9 |

### reading_time_continuous

A continuous predictor with a lognormal reading-time outcome. The
`lognormal` family, 4000 rows.

| subject | item | narration |  SyntaxPC | CoherencePC |       age | reading_time_per_word |
|--------:|-----:|:----------|----------:|------------:|----------:|----------------------:|
|       1 |    1 | off       | -1.739797 |  -0.0215117 | 0.0632888 |               0.37485 |
|       1 |    1 | on        | -1.739797 |  -0.0215117 | 0.0632888 |               0.32765 |
|       1 |    2 | off       | -1.891749 |   0.7268579 | 0.0632888 |               0.15592 |
|       1 |    2 | on        | -1.891749 |   0.7268579 | 0.0632888 |               0.23166 |
