# A single sweep over any field of a design specification.
#
# power_curve_mixed() and precision_curve() were the same six lines: copy the specification, set
# units$subject$n, call the analysis, bind the rows. Sample size was therefore the only axis a user
# could vary without writing the loop themselves, and effect size, which is the axis a design
# analysis usually needs next, was not reachable. One sweep over an addressed field covers both, and
# the two curve functions become thin wrappers that keep working.

# Set a nested field, addressed either as "units$subject$n" or as c("units", "subject", "n").
.set_path <- function(spec, path, value) {
  keys <- if (length(path) == 1L) strsplit(path, "$", fixed = TRUE)[[1]] else as.character(path)
  keys <- keys[nzchar(keys)]
  if (!length(keys)) stop("`path` must name at least one field", call. = FALSE)
  .assign <- function(x, k, v) {
    if (length(k) == 1L) { x[[k[1]]] <- v; return(x) }
    if (is.null(x[[k[1]]]))
      stop("`path` names '", k[1], "', which this specification does not have", call. = FALSE)
    x[[k[1]]] <- .assign(x[[k[1]]], k[-1], v)
    x
  }
  .assign(spec, keys, value)
}

# Read a nested field: NULL when a path does not exist, so that sweep_spec() can say so, and
# otherwise the value a design_conditions() condition is merged into.
.get_path <- function(spec, path) {
  keys <- if (length(path) == 1L) strsplit(path, "$", fixed = TRUE)[[1]] else as.character(path)
  keys <- keys[nzchar(keys)]
  cur <- spec
  for (k in keys) {
    if (!is.list(cur) || is.null(cur[[k]])) return(NULL)
    cur <- cur[[k]]
  }
  cur
}

#' Sweep an analysis over one field of a design specification
#'
#' Vary a single field of a specification across a set of values, run an analysis at each, and bind
#' the results into one data frame. Sample size is the axis users sweep most often, and
#' [power_curve_mixed()] and [precision_curve()] are wrappers around this for it, but any field can
#' be swept, including an effect size, a random-effect standard deviation, a residual standard
#' deviation, or the number of items per subject.
#'
#' @details
#' The specification is validated once, before the sweep, so a mistake in it is reported before any
#' fitting starts rather than repeated at every grid point.
#'
#' A value may be a scalar, which replaces the addressed field, or a list, which replaces it
#' wholesale. The conditions [design_conditions()] builds are the exception, because each names
#' only the effects it varies. Each condition is merged into the addressed list with
#' [utils::modifyList()], so the coefficients it names take its values and every other
#' coefficient keeps the specification's, in the specification's order. That is the order in
#' which the linear predictor is summed, so a condition simulates exactly the data of the
#' specification edited by hand. A condition may name only coefficients the specification
#' already has, and any other name is refused.
#'
#' The result of `fn` is coerced to a data frame: a `pilotr_power` object becomes one row per focal
#' effect, a data frame is used as it stands, and a plain list becomes a single row. The swept value
#' is added as a leading column named by `.name`. A list has no single value to record, so a sweep
#' over lists, the conditions [design_conditions()] builds included, records the grid index there.
#' [solve_curve()] then needs `x` to name the column that holds the effect. To vary a single
#' effect, sweep its own coefficient. For example,
#' `sweep_spec(spec, "fixed$coefficients$cond", c(0, 0.02, 0.04), power_mixed)` keeps the other
#' coefficients and names the leading column `cond`, so
#' `solve_curve(curve, target = 0.8, effect = "cond", transform = "identity")` returns the
#' minimum detectable effect.
#'
#' @param spec A design specification (path or list).
#' @param path The field to vary, as `"units$subject$n"` or `c("units", "subject", "n")`.
#' @param values A vector or list of values to set the field to, one grid point each.
#'   [design_conditions()] builds such a list for an effect-size sweep.
#' @param fn The analysis to run at each grid point, for example [power_mixed()] or
#'   [precision_design()]. It is called as `fn(spec, ...)`.
#' @param ... Further arguments passed to `fn` at every grid point.
#' @param .name Name for the column recording the swept value. Defaults to the last element of
#'   `path`.
#' @return A data frame binding the results, with the swept value as the leading column, or the
#'   grid index when the values are lists. [solve_curve()] reads that leading column by default,
#'   so a sweep over single values goes straight into a solve.
#' @examples
#' \donttest{
#' if (requireNamespace("lme4", quietly = TRUE) &&
#'     requireNamespace("lmerTest", quietly = TRUE)) {
#'   spec <- build_spec(list(name = "p", seed = 1, design_kind = "within",
#'     include_items = TRUE, n_subject = 12, n_item = 12, factor_name = "cond",
#'     lev1 = "a", lev2 = "b", intercept = 6, effect = 0.05,
#'     subj_int_sd = 0.12, subj_slope_sd = 0.04, subj_corr = 0.2,
#'     item_int_sd = 0.08, item_slope_sd = 0.02, item_corr = -0.1,
#'     family = "shifted_lognormal", resp_name = "", sigma = 0.3, shift = 200))
#'
#'   # Sample size, the same sweep power_curve_mixed() performs.
#'   sweep_spec(spec, "units$subject$n", c(12, 18), power_mixed, n_sims = 8)
#'
#'   # Effect size, which the old curve functions could not reach.
#'   sweep_spec(spec, "fixed$coefficients",
#'              design_conditions(effect = c(0.03, 0.06)), power_mixed, n_sims = 8)
#' }
#' }
#' @seealso [design_conditions()] to build effect-size grids, [power_mixed()] and
#'   [precision_design()] for the analyses usually swept, and [solve_curve()] to solve the
#'   resulting curve for the swept value that meets a target.
#' @export
sweep_spec <- function(spec, path, values, fn, ..., .name = NULL) {
  spec <- .as_spec(spec)
  if (!is.function(fn)) stop("`fn` must be a function", call. = FALSE)
  if (!length(values)) stop("`values` must contain at least one value", call. = FALSE)
  current <- .get_path(spec, path)
  if (is.null(current))
    stop("`path` does not address a field of this specification: ",
         paste(path, collapse = "$"), call. = FALSE)
  merge <- inherits(values, "pilotr_conditions")
  if (merge) .check_conditions(values, current, path)
  if (is.null(.name)) {
    keys <- if (length(path) == 1L) strsplit(path, "$", fixed = TRUE)[[1]] else as.character(path)
    .name <- keys[length(keys)]
  }

  scalar_values <- !is.list(values)
  parts <- lapply(seq_along(values), function(i) {
    value <- if (is.list(values)) values[[i]] else values[i]
    # A condition names only the effects it varies, so it is merged into the field. modifyList()
    # keeps the field's own order, which for the coefficients is the order in which the linear
    # predictor is summed. Reversing it moves 691 of the 4,000 unrounded responses of
    # reading_time_continuous by up to 3.3e-16.
    if (merge) value <- utils::modifyList(current, value)
    # The specification is already validated, and each grid point only replaces one field, so
    # revalidating at every point would repeat the same work; fn() validates its own argument.
    res <- fn(.set_path(spec, path, value), ...)
    frame <- if (inherits(res, "pilotr_power")) .as_power_frame(res)
      else if (is.data.frame(res)) res
      else as.data.frame(res, stringsAsFactors = FALSE)
    lead <- if (scalar_values) value else i
    cbind(stats::setNames(data.frame(lead), .name), frame)
  })
  do.call(rbind, parts)
}

# Conditions may name only what the addressed list already holds. Merging another name would
# append a coefficient the specification never declared, at the end of its order, and the
# validator accepts any such key that resolves to a column.
.check_conditions <- function(values, current, path) {
  if (!is.list(current) || (length(current) && is.null(names(current))))
    stop("design_conditions() merges its conditions into a named list such as ",
         "fixed$coefficients, and `path` addresses ", paste(path, collapse = "$"),
         ", which is not one", call. = FALSE)
  have <- names(current)
  unknown <- setdiff(unique(unlist(lapply(values, names), use.names = FALSE)), have)
  if (length(unknown))
    stop(sprintf(
      "design_conditions() names %s, which %s not %s of this specification; its coefficients are %s",
      paste(sprintf("'%s'", unknown), collapse = ", "),
      if (length(unknown) > 1L) "are" else "is",
      if (length(unknown) > 1L) "coefficients" else "a coefficient",
      if (length(have)) paste(sprintf("'%s'", have), collapse = ", ") else "(none)"),
      call. = FALSE)
  invisible(NULL)
}

#' Build a grid of fixed-effect coefficient sets
#'
#' Produce the conditions needed to sweep an effect size with [sweep_spec()], including a
#' condition in which every named effect is zero, so that the same run shows both what a design
#' detects and how often it declares something when there is nothing to find.
#'
#' @details
#' Named arguments give the values each effect should take, and are recycled to a common length, so
#' `design_conditions(cond = c(0.02, 0.05), age = 0.1)` produces two conditions, both with `age` at
#' 0.1. The all-zero condition comes first and is shared, since the Type I error rate is a property
#' of the design rather than of any one effect size.
#'
#' [sweep_spec()] merges each condition into the specification's own coefficients and keeps their
#' order, so any coefficient not named here keeps its value and a sweep varies only the effects it
#' names. It refuses a condition that names a coefficient the specification lacks. Subsetting the
#' result with `[` keeps its class, and so does `c()` when the result comes first. The three
#' conditions of `c(design_conditions(cond = 0.02), design_conditions(age = 0.1, .null = FALSE))`
#' are therefore merged in the same way. A list assembled by other means replaces the
#' coefficients wholesale.
#'
#' @param ... Named numeric vectors, one per coefficient to vary.
#' @param .null Whether to prepend a condition with every named effect set to zero. `TRUE` by
#'   default.
#' @param .base Optional named list of coefficients to include in every condition, the null
#'   condition among them. [sweep_spec()] already keeps the coefficients a condition does not
#'   name. `.base` is therefore needed only to hold a coefficient at a value other than the
#'   specification's, or for a condition used outside [sweep_spec()], such as one assigned to
#'   `fixed$coefficients` by hand.
#' @return A list of class `pilotr_conditions` holding one named list of coefficients per
#'   condition. The class is how [sweep_spec()] recognises conditions to merge.
#' @examples
#' design_conditions(effect = c(0.03, 0.06))
#' design_conditions(cond = c(0.02, 0.05), age = 0.1, .null = FALSE)
#' @seealso [sweep_spec()], which consumes this.
#' @export
design_conditions <- function(..., .null = TRUE, .base = NULL) {
  effects <- list(...)
  if (!length(effects) || is.null(names(effects)) || any(!nzchar(names(effects))))
    stop("name every effect, as in design_conditions(cond = c(0.02, 0.05))", call. = FALSE)
  if (!all(vapply(effects, is.numeric, logical(1))))
    stop("every effect must be numeric", call. = FALSE)
  n <- max(vapply(effects, length, integer(1)))
  effects <- lapply(effects, function(v) rep(v, length.out = n))

  conditions <- lapply(seq_len(n), function(i) {
    cond <- lapply(effects, function(v) v[i])
    utils::modifyList(as.list(.base %||NULL% list()), cond)
  })
  if (isTRUE(.null)) {
    zero <- lapply(effects, function(v) 0)
    conditions <- c(list(utils::modifyList(as.list(.base %||NULL% list()), zero)), conditions)
  }
  structure(conditions, class = c("pilotr_conditions", "list"))
}

# Without this method `[` returns a plain list, which sweep_spec() substitutes wholesale, so a
# subset such as conditions[-1] would again drop every coefficient it does not name. Base R's
# `[.simple.list` keeps its class in the same way.
#' @export
#' @noRd
`[.pilotr_conditions` <- function(x, i, ...) structure(NextMethod("["), class = class(x))

# Joining two sets of conditions with c(), to vary one effect and then another in the same run,
# would likewise return a plain list. Base R's `c.noquote` keeps its class in the same way. c()
# dispatches on its first argument, so the conditions have to come first.
#' @export
#' @noRd
c.pilotr_conditions <- function(..., recursive = FALSE) {
  out <- NextMethod("c")
  if (is.list(out)) class(out) <- c("pilotr_conditions", "list")
  out
}

# The conditions print as the plain list they hold, without the class attribute that
# sweep_spec() reads.
#' @export
#' @noRd
print.pilotr_conditions <- function(x, ...) {
  print(unclass(x), ...)
  invisible(x)
}

# A local null-coalesce, named so that it cannot be confused with base R's `%||%` (added in 4.4)
# and so that pilotr keeps working on older R.
`%||NULL%` <- function(a, b) if (is.null(a)) b else a
