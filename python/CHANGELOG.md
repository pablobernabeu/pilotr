# Changelog

All notable changes to this project are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2026-07-09

### Added

- `power`, `power_mixed` and `power_curve` gain a `workers` argument that spreads the
  Monte Carlo replicates across local processes with
  `concurrent.futures.ProcessPoolExecutor`. Because every replicate seeds the shared
  cross-language RNG from its own index, any worker count returns results identical to
  a serial run.
- `power_curve` starts one process pool and reuses it across the whole sweep, so a
  high-resolution curve pays the pool start-up cost only once.

## [0.1.0] - 2026-07-08

The first public version of the Python package. It simulates experimental and
behavioural data from a portable JSON design specification, producing output that
is bit-identical to the R package of the same name given the same specification
and seed.

### Added

- `simulate` generates a dataset from a design specification supplied as a plain
  dictionary or a JSON file, with `load_spec` reading a specification authored
  elsewhere, such as one downloaded from the no-code app.
- Response families for Gaussian outcomes, reaction times, accuracy, counts,
  ordinal responses and Beta-distributed proportions.
- Power and design analysis through `power`, `power_curve` and `power_mixed`,
  with `scipy` and `statsmodels` supplied as optional extras so the generative
  core stays dependency-free.
- A cross-language random-number contract that reproduces the R package's data to
  full floating-point precision.
- A documentation site whose guides execute their examples at build time, so the
  tables and figures shown are real `pilotr` output.
