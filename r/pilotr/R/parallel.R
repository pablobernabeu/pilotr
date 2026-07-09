# Internal helpers for parallelising the Monte Carlo replicate loops. Every replicate
# seeds the shared cross-language RNG from its own index (base seed + index), so the
# results are identical whether the replicates run serially or on a PSOCK cluster, and
# whatever the worker count. The helpers therefore only distribute work; they never
# touch the random-number state.

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
