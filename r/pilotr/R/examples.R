# Locate the bundled example specifications. The JSON files under inst/examples/
# are copies of the repository's spec/examples/, shipped inside the package so
# they are reachable with system.file() from an installed copy. They are the
# same files that drive the Python twin and the no-code app, so a design authored
# once runs unchanged across all three.

#' Example design specifications shipped with pilotr
#'
#' pilotr ships one ready-to-run specification per design family, as JSON, in the
#' package's `examples/` directory. These are the same files that drive the
#' Python twin and the no-code app, so a design authored in one place runs
#' unchanged in the others. `pilotr_example()` lists them, or returns the path to
#' one for [load_spec()].
#'
#' @param name The base name of an example, with or without the `.json`
#'   extension, for example `"between_2group_gaussian"`. When `NULL` (the
#'   default), the available example names are returned instead of a path.
#' @return When `name` is `NULL`, a character vector of the available example
#'   names. Otherwise, the full path to that example's JSON file, ready to pass
#'   to [load_spec()].
#' @seealso [load_spec()] to read a specification and
#'   [simulate_design()] to simulate from it.
#' @examples
#' pilotr_example()                       # the available examples
#' spec <- load_spec(pilotr_example("between_2group_gaussian"))
#' head(simulate_design(spec))
#' @export
pilotr_example <- function(name = NULL) {
  dir <- system.file("examples", package = "pilotr")
  available <- sort(tools::file_path_sans_ext(
    list.files(dir, pattern = "[.]json$")
  ))
  if (is.null(name)) {
    return(available)
  }
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("`name` must be a single example name, or NULL to list them.",
         call. = FALSE)
  }
  base <- tools::file_path_sans_ext(name)
  if (!base %in% available) {
    stop(sprintf("Unknown example '%s'. Available: %s.",
                 name, paste(available, collapse = ", ")), call. = FALSE)
  }
  system.file("examples", paste0(base, ".json"), package = "pilotr")
}
