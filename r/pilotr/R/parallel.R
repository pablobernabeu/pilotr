# Internal helpers for parallelising the Monte Carlo replicate loops. Every replicate takes its
# own seed from replicate_seeds(), so the results are identical whether the replicates run
# serially or on a PSOCK cluster, and whatever the worker count. The helpers therefore only
# distribute work; they never touch the random-number state.

# Validate the `workers` argument: a single positive whole number.
.check_workers <- function(workers) {
  if (!is.numeric(workers) || length(workers) != 1L || is.na(workers) ||
      workers < 1 || workers != round(workers))
    stop("`workers` must be a single positive whole number.", call. = FALSE)
  as.integer(workers)
}

# Apply FUN over the replicate indices: serially when `cl` is NULL, otherwise with
# load balancing over the cluster. parLapplyLB returns results in input order, so the
# downstream reductions match the serial code exactly.
.p_lapply <- function(indices, FUN, cl = NULL, ...) {
  if (is.null(cl)) lapply(indices, FUN, ...)
  else parallel::parLapplyLB(cl, indices, FUN, ...)
}

#' Seeds for the replicates of a Monte Carlo run
#'
#' The seeds pilotr's power and precision loops give to their replicates, derived from a
#' specification's own seed. Exported so that a hand-written replicate loop, or a cluster array
#' task that has to reproduce one replicate on its own, can use the same rule and so land on the
#' same data.
#'
#' @details
#' Until 0.3 the rule was `base + (i - 1)`. Consecutive seeds are not independent streams in this
#' generator: seeding sets `s1` to `1 + (seed mod 2147483562)` and `s2` from `s1`, and only ten
#' warm-up draws are discarded, so replicate `i` and replicate `i + 1` begin a few steps apart in
#' the same sequence rather than in unrelated parts of it. Measured over 2,000 replicates, the
#' first draw of replicate `i` correlated 0.95 with the first draw of replicate `i + 1`.
#'
#' An arithmetic scramble does not fix that. Adding a Weyl increment and applying a Lehmer step
#' leaves the seeds in arithmetic progression with a longer stride, and since the seeding rule is
#' itself linear in the seed, the first draws remained correlated at -0.27. The problem is
#' linearity, so no linear remedy addresses it.
#'
#' Drawing the seeds from the shared generator does work. Successive outputs of the combined
#' generator are what that generator exists to make look independent, so the seeds inherit it: the
#' same measurement gives -0.02, and a Ljung-Box test over the resulting replicate means moves from
#' p below 0.0001 to p of 0.94. The whole vector costs `n` draws computed once rather than any work
#' per replicate. Duplicates are skipped, so no two replicates are handed the same seed and
#' silently produce identical data, and the skipping is deterministic, so the R and 'Python'
#' implementations still agree.
#'
#' This changed every number pilotr produced before 0.3, and the first replicate no longer uses the
#' specification's own seed. Both are deliberate; pin an earlier version to reproduce earlier
#' output.
#'
#' @param base The specification's seed.
#' @param n How many replicate seeds to return.
#' @return A numeric vector of `n` distinct seeds.
#' @examples
#' replicate_seeds(90210, 5)
#'
#' # The same data as replicate 3 of a power run over this specification.
#' spec <- build_spec(list(name = "d", seed = 90210, design_kind = "between",
#'   factor_name = "g", lev1 = "a", lev2 = "b", n_subject = 20,
#'   intercept = 0, effect = 0.5, family = "gaussian", resp_name = "", sigma = 1))
#' spec$seed <- replicate_seeds(90210, 3)[3]
#' head(simulate_design(spec), 3)
#' @export
replicate_seeds <- function(base, n) {
  rng <- make_rng(base)
  out <- numeric(n)
  seen <- new.env(hash = TRUE, parent = emptyenv())
  k <- 0L
  while (k < n) {
    s <- floor(rng$uniform() * (.M1 - 1)) + 1
    key <- format(s, scientific = FALSE)
    if (!exists(key, envir = seen, inherits = FALSE)) {
      k <- k + 1L
      out[k] <- s
      assign(key, TRUE, envir = seen)
    }
  }
  out
}
