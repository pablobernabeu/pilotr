# Load the pilotr engine from source, for the scripts under examples/ that
# exercise it without installing the package.
#
# It loads every file in R/, so that no script has to name the engine files it
# needs and keep that list in step with R/ by hand. A list that leaves a file out
# goes unnoticed until a function defined there is called.
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
# examples/ is one level below the package root.
pilotr_engine_dir <- function(script_dir) file.path(script_dir, "..", "R")
