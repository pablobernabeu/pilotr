# Contributing to pilotr

Thank you for considering a contribution. Issues and pull requests are both
welcome, whether they fix a bug, improve the documentation or add a feature.

## Reporting a problem or suggesting a feature

Please open an issue at <https://github.com/pablobernabeu/pilotr/issues>. A small
reproducible example helps a great deal, and for a bug it is useful to include a
minimal design specification and the output of `sessionInfo()`.

## Setting up for development

Clone the repository and install the development dependencies:

```r
# install.packages("pak")
pak::pak("pablobernabeu/pilotr")
pak::pak(c("devtools", "roxygen2", "testthat", "spelling"))
```

The package is developed with the usual devtools workflow:

```r
devtools::document()   # regenerate man/ and NAMESPACE after editing roxygen
devtools::test()       # run the test suite
devtools::check()      # a full R CMD check
```

## How the code is organised

pilotr and its Python twin read the same JSON design specification and share a
cross-language random-number generator, so that a given specification and seed
produce bit-identical data in both languages. A change to the simulation core
should keep that guarantee, and a pull request that touches it is easiest to
accept when it adds a test showing the two implementations still agree.

The tests run entirely offline. The prose follows British spelling, and running
`spelling::spell_check_package()` before a pull request keeps the word list tidy.

## Submitting a pull request

Please base your work on `main`, keep the change focused, and add or update tests
and documentation alongside the code. Running `devtools::document()`,
`devtools::test()` and `devtools::check()` locally before opening the pull
request saves a round trip.

By contributing you agree that your contribution is licensed under the same MIT
licence as the package.
