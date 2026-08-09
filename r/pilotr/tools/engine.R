# Load the pilotr engine from source, for the scripts under examples/ and tests/
# that exercise it without installing the package.
#
# Those scripts each used to name the engine files they needed. The lists had to
# be kept in step with R/ by hand, once per script, and they were not: the 0.3.0
# work added validate.R, which load_spec() and build_spec() both reach into, and
# every list that omitted it kept working until the moment it called one of them,
# then failed with "could not find function \"validate_spec\"". Loading the
# directory removes the lists, and with them the thing that drifts.
#
# Every file in R/ defines functions and constants and nothing else -- no
# library() call, no other top-level side effect -- so loading them all is cheap
# and the order does not matter. pilotr-package.R holds only the roxygen
# "_PACKAGE" sentinel and has nothing to define.

load_pilotr_engine <- function(r_dir) {
  if (!dir.exists(r_dir)) {
    stop("No engine directory at ", r_dir, ".", call. = FALSE)
  }
  files <- sort(list.files(r_dir, pattern = "[.][Rr]$", full.names = TRUE))
  files <- files[basename(files) != "pilotr-package.R"]
  if (!length(files)) {
    stop("No engine sources found in ", r_dir, ".", call. = FALSE)
  }
  # source() into the global environment, which is where these scripts expect
  # the engine to land; the default local = FALSE does that from anywhere.
  for (f in files) source(f)
  invisible(normalizePath(files))
}

# Resolve the engine directory from the directory holding the calling script.
# examples/ and tests/ are both one level below the package root, so the path is
# the same from either.
pilotr_engine_dir <- function(script_dir) file.path(script_dir, "..", "R")
