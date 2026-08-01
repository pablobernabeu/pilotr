# Cross-language parity harness, R side.
#
# Simulates every specification in spec/examples/ from the *local* package sources (not an
# installed copy) and writes one canonical dump per case. The Python driver writes the same
# dumps from python/pilotr/, and compare.py diffs them byte for byte.
#
# Each specification is run twice: "asis", and "noround" with `response.round` deleted. The
# second case matters because six of the eight shipped examples set `round`, which quantises
# away exactly the last-ulp divergences this harness exists to detect.
#
# Usage: Rscript tools/parity/run_r.R [output-dir]

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "..", ".."), mustWork = FALSE)
if (!nzchar(root) || is.na(root)) root <- normalizePath(".")
out_dir <- if (length(args) >= 1) args[1] else file.path(root, "tools", "parity", "out", "r")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Source the whole package rather than a hand-listed subset, so that a new file cannot leave the
# harness silently running against a stale definition. The generative path pulls in the
# validator, so every case is validated here too.
src <- file.path(root, "r", "pilotr", "R")
for (f in sort(list.files(src, pattern = "\\.R$", full.names = TRUE))) source(f)

# printf must round-trip a double exactly, or the whole comparison is meaningless. Check the
# formatter itself before trusting it on simulated data.
.check_formatter <- function() {
  probes <- c(1/3, 0.1, pi, 1e-300, .Machine$double.xmax, -2.5e-17, 0, 1, -1,
              .Machine$double.eps, 1 - .Machine$double.eps / 2)
  for (v in probes) {
    if (!identical(as.numeric(sprintf("%.17g", v)), v))
      stop("sprintf('%.17g') does not round-trip ", format(v, digits = 22), call. = FALSE)
  }
}
.check_formatter()

# Language-neutral cell rendering: character stays verbatim, everything else becomes a double
# rendered at 17 significant digits. Integers therefore print as "1", matching Python.
.cell <- function(v) if (is.character(v) || is.factor(v)) as.character(v) else sprintf("%.17g", as.numeric(v))

.dump <- function(d, path) {
  con <- file(path, open = "wb")           # binary: LF endings on every platform
  on.exit(close(con))
  lines <- character(nrow(d) + 1L)
  lines[1] <- paste(names(d), collapse = ",")
  cols <- lapply(d, function(col) vapply(col, .cell, character(1), USE.NAMES = FALSE))
  for (i in seq_len(nrow(d)))
    lines[i + 1L] <- paste(vapply(cols, `[`, character(1), i), collapse = ",")
  writeLines(lines, con, sep = "\n")
}

# Two sources. spec/examples/ is what ships, and must never change silently.
# tools/parity/cases/ holds adversarial specifications that deliberately exercise what the
# shipped examples do not: long random-effect vectors (so the Cholesky and matrix-vector
# inner products run to three terms or more, where the two languages' built-in reductions
# diverge), long linear predictors, interaction random slopes, and a Gaussian response, so
# that no libm transcendental masks or manufactures a difference.
spec_dirs <- c(file.path(root, "spec", "examples"), file.path(root, "tools", "parity", "cases"))
specs <- unlist(lapply(spec_dirs, function(d)
  if (dir.exists(d)) sort(list.files(d, pattern = "\\.json$", full.names = TRUE)) else character(0)))
if (!length(specs)) stop("no specifications found in ", paste(spec_dirs, collapse = " or "), call. = FALSE)

for (p in specs) {
  base <- sub("\\.json$", "", basename(p))
  for (variant in c("asis", "noround")) {
    spec <- load_spec(p)
    if (variant == "noround") spec$response$round <- NULL
    d <- simulate_design(spec)
    .dump(d, file.path(out_dir, paste0(base, ".", variant, ".txt")))
  }
  cat("simulated:", base, "\n")
}

cat("R dumps written to", out_dir, "\n")
