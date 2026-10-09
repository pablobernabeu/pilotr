# Validate a design specification

Check a design specification against the portable schema and against the
cross-field rules the schema cannot express, and check that its declared
`spec_version` is one this implementation understands. Called by
[`load_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/load_spec.md)
by default.

## Usage

``` r
validate_spec(spec, strict = TRUE)
```

## Arguments

- spec:

  A design specification (path or list).

- strict:

  Whether an unrecognised field is an error (the default) or a warning.
  Set `FALSE` to load a specification carrying private annotations,
  accepting that a misspelled field will then be ignored in silence.

## Value

The specification, invisibly, so that the call can be chained.

## Details

Validation exists because several ways of getting a specification wrong
produce plausible data and no error at all. A mistyped coefficient key
resolves to no column and so silently sets that effect to zero, which
generates exactly the data of a null design and reports success. A
response parameter left over from another family is ignored. Neither is
detectable in the output, which is why they are refused here.

Names are checked together as well as one by one. Two columns with one
name, a contrast column defined by two factors or named like another
column, an interaction whose analysis column already exists, a level
listed twice, a correlation that pairs a term with itself or gives one
pair twice, a factor both between and within a unit, and a `random.item`
entry in a design without items each validated and then changed the data
in silence. The rules are set out under Names in the specification. A
within factor whose `vary_within` leaves out a unit of the design draws
a warning in either mode, since pilotr crosses it with every unit
anyway.

At most one factor may be between each unit. pilotr assigns the levels
of each between factor to blocks of its unit on its own, so two between
factors over one unit fall into the same or overlapping blocks. A 2 x 2
between-subjects design over 40 subjects gave cells of 20, 0, 0 and 20,
and its effects could not be estimated apart. Such a design is written
as one between factor whose levels are the cells, as the specification
shows under Worked encodings, and
[`spec_from_model()`](https://pablobernabeu.github.io/pilotr/r/reference/spec_from_model.md)
builds that factor from a fitted pilot.

Version negotiation covers the other direction. A specification that
uses a feature introduced in 0.3 is read differently by a 0.2
implementation, so it must declare 0.3 or later. A specification
declaring a version newer than this implementation is refused outright.
A specification with no `spec_version` is treated as 0.2, which is what
every specification written before the field existed is.

Every field is read by its exact name, as the 'Python' twin reads it,
and a list that repeats a name within one object is refused, as
[`load_spec()`](https://pablobernabeu.github.io/pilotr/r/reference/load_spec.md)
refuses a file that does. A seed must be a whole number within plus or
minus 2^53 - 1, the range in which both twins read a JSON integer
exactly.

## Examples

``` r
spec <- build_spec(list(name = "demo", seed = 1, design_kind = "between",
  factor_name = "group", lev1 = "a", lev2 = "b", n_subject = 20,
  intercept = 0, effect = 0.5, family = "gaussian", resp_name = "", sigma = 1))
validate_spec(spec)

# A mistyped coefficient key is refused, where it used to pass as a zero effect.
bad <- spec
bad$fixed$coefficients <- list(effct = 0.5)
try(validate_spec(bad))
#> Error : invalid design specification:
#>   - fixed.coefficients 'effct' names 'effct', which is neither a contrast column nor a predictor; available columns are 'effect'. An unresolved key contributes zero, so this would silently drop the term
```
