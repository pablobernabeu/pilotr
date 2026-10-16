# Seeds for the replicates of a Monte Carlo run

The seeds pilotr's power and precision loops give to their replicates,
derived from a specification's own seed. Exported so that a hand-written
replicate loop, or a cluster array task that has to reproduce one
replicate on its own, can use the same rule and so land on the same
data.

## Usage

``` r
replicate_seeds(base, n)
```

## Arguments

- base:

  The specification's seed.

- n:

  How many replicate seeds to return.

## Value

A numeric vector of `n` distinct seeds.

## Details

Until 0.3 the rule was `base + (i - 1)`, and consecutive seeds are not
independent streams in this generator. Seeding is affine in the seed,
setting `s1` to `1 + (seed mod 2147483562)` and `s2` from `s1`, and each
component recurrence is linear, so at every position `k` the draws of
consecutive seeds differ by very nearly the same amount mod 1. For the
first draw that step is -0.0091: seeds 1000 to 1005 begin with 0.7027,
0.6937, 0.6846, 0.6755, 0.6665 and 0.6574. Draws in arithmetic
progression mod 1 with step `d` correlate `1 - 6|d|(1 - |d|)`, so over
2,000 replicates the first draw of replicate `i` correlated 0.95 with
the first draw of replicate `i + 1`. A longer warm-up changes `d` and
leaves the progression in place: with 50 warm-up draws in place of 10,
the step is -0.32 and the correlation -0.30.

An arithmetic scramble does not fix that. Adding a Weyl increment and
applying a Lehmer step leaves the seeds in arithmetic progression with a
longer stride, and since the seeding rule is affine in the seed, the
first draws remained correlated at -0.27. Any affine rule for the seeds
keeps the draws in progression, so no linear remedy addresses it.
Matsumoto et al. (2007) describe the same defect in the initialisation
of other linear generators.

Drawing the seeds from the shared generator does work. Successive
outputs of the combined generator are what that generator exists to make
look independent, so the seeds inherit it: the same measurement gives
-0.02, and a Ljung-Box test over the resulting replicate means moves
from p below 0.0001 to p of 0.94. The whole vector costs `n` draws,
computed once for the loop. Duplicates are skipped, so no two replicates
are handed the same seed and silently produce identical data, and the
skipping is deterministic, so the R and 'Python' implementations still
agree.

This changed every number pilotr produced before 0.3, and the first
replicate no longer uses the specification's own seed. Both are
deliberate; pin an earlier version to reproduce earlier output.

## References

Matsumoto, M., Wada, I., Kuramoto, A. and Ashihara, H. (2007). Common
defects in initialization of pseudorandom number generators. *ACM
Transactions on Modeling and Computer Simulation*, 17(4), 15.
[doi:10.1145/1276927.1276928](https://doi.org/10.1145/1276927.1276928)

## Examples

``` r
replicate_seeds(90210, 5)
#> [1] 1222973854 1397434102  271722066 2119005144 1916898993

# The same data as replicate 3 of a power run over this specification.
spec <- build_spec(list(name = "d", seed = 90210, design_kind = "between",
  factor_name = "g", lev1 = "a", lev2 = "b", n_subject = 20,
  intercept = 0, effect = 0.5, family = "gaussian", resp_name = "", sigma = 1))
spec$seed <- replicate_seeds(90210, 3)[3]
head(simulate_design(spec), 3)
#>   subject g   score
#> 1       1 a -1.1770
#> 2       2 a  0.6361
#> 3       3 a  0.8406
```
