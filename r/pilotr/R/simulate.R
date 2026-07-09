# Specification parsing and the generative engine. This is a bit-identical mirror of
# pilotr/simulate.py. It supports categorical contrasts and continuous predictors (with
# interactions) as fixed effects, together with crossed by-subject and by-item random
# intercepts and slopes (on either unit).

#' Load a design specification from a JSON file
#'
#' @param path Path to a JSON design-specification file.
#' @return The specification as a nested list, with sub-lists left unsimplified so that the
#'   structure round-trips exactly. Pass the result to [simulate_design()].
#' @examples
#' spec <- build_spec(list(name = "demo", seed = 1, design_kind = "between",
#'   factor_name = "group", lev1 = "a", lev2 = "b", n_subject = 20,
#'   intercept = 0, effect = 0.5, family = "gaussian", resp_name = "", sigma = 1))
#' f <- tempfile(fileext = ".json")
#' writeLines(spec_json(spec), f)
#' identical(simulate_design(load_spec(f)), simulate_design(spec))
#' @export
load_spec <- function(path) {
  jsonlite::fromJSON(path, simplifyVector = TRUE, simplifyDataFrame = FALSE,
                     simplifyMatrix = FALSE)
}

# Follows itertools.product order, in which the last index varies fastest. Returns a list
# of 0-based index vectors.
.product_indices <- function(sizes) {
  if (length(sizes) == 0) return(list(integer(0)))
  n <- length(sizes); idx <- rep(0L, n); out <- list()
  repeat {
    out[[length(out) + 1L]] <- idx
    k <- n
    while (k >= 1) {
      idx[k] <- idx[k] + 1L
      if (idx[k] < sizes[k]) break
      idx[k] <- 0L; k <- k - 1L
    }
    if (k < 1) break
  }
  out
}

# Returns the columns and the lower-Cholesky factor L for one unit's random-effect covariance.
.ranef <- function(u) {
  slope_names <- names(u$slopes)
  cols <- c("intercept", slope_names)
  sds <- c(u$intercept_sd, vapply(slope_names, function(s) u$slopes[[s]], numeric(1)))
  n <- length(cols)
  R <- diag(n)
  if (!is.null(u$correlations)) for (key in names(u$correlations)) {
    parts <- trimws(strsplit(gsub("~", ",", key), ",")[[1]])
    i <- match(parts[1], cols); j <- match(parts[2], cols)
    R[i, j] <- R[j, i] <- u$correlations[[key]]
  }
  cov <- outer(sds, sds) * R
  list(cols = cols, L = .cholesky(cov))
}

# A coefficient/slope key is a column name or an 'a:b' interaction (product of columns).
.design_value <- function(cvals, key) {
  if (grepl(":", key, fixed = TRUE)) {
    v <- 1
    for (pp in strsplit(key, ":", fixed = TRUE)[[1]]) v <- v * (if (is.null(cvals[[pp]])) 0 else cvals[[pp]])
    return(v)
  }
  if (is.null(cvals[[key]])) 0 else cvals[[key]]
}

# Samples m distinct items from 1..n_items via a partial Fisher-Yates shuffle on the shared
# RNG, which produces partial crossing. This is bit-identical with simulate.py's _sample_items.
.sample_items <- function(rng, n_items, m) {
  pool <- 1:n_items
  for (k in 0:(m - 1)) {
    j <- k + floor(rng$uniform() * (n_items - k))
    tmp <- pool[k + 1]; pool[k + 1] <- pool[j + 1]; pool[j + 1] <- tmp
  }
  sort(pool[1:m])
}

#' Simulate a data set from a design specification
#'
#' Generate an analysis-ready data set from a portable design specification: a linear
#' predictor built from fixed effect sizes (categorical contrasts, continuous predictors,
#' and their interactions) plus crossed by-subject and by-item random intercepts and slopes,
#' mapped through the chosen response family.
#'
#' @param spec A design specification, given either as a path to a JSON file or as an
#'   already-parsed list (for example from [build_spec()] or [load_spec()]).
#' @return A data frame with one row per observation, containing a `subject` column, an
#'   optional `item` column, any grouping, factor, and continuous-predictor columns, and the
#'   response column named by the specification.
#' @examples
#' spec <- build_spec(list(name = "demo", seed = 1, design_kind = "between",
#'   factor_name = "group", lev1 = "control", lev2 = "treatment", n_subject = 40,
#'   intercept = 100, effect = 5, family = "gaussian", resp_name = "", sigma = 10))
#' head(simulate_design(spec))
#' @export
simulate_design <- function(spec) {
  if (is.character(spec)) spec <- load_spec(spec)

  S <- spec$units$subject$n
  has_item <- !is.null(spec$units$item)
  I <- if (has_item) spec$units$item$n else 1L

  factors <- spec$factors
  predictors <- spec$predictors
  within <- Filter(function(f) !is.null(f$vary_within), factors)
  between <- Filter(function(f) !is.null(f$between), factors)
  within_sizes <- vapply(within, function(f) length(f$levels), integer(1))

  rng <- make_rng(spec$seed)
  per_subject <- if (has_item) spec$units$item$per_subject else NULL
  if (!is.null(per_subject)) {
    if (per_subject < 1)
      stop("per_subject (", per_subject, ") must be at least 1", call. = FALSE)
    if (per_subject > I)
      stop("per_subject (", per_subject, ") cannot exceed the number of items (", I, ")",
           call. = FALSE)
  }

  # ---- canonical row order; per-subject item subsets (if any) are the first RNG draws ----
  rows <- list()
  for (s in 1:S) {
    items_s <- if (!is.null(per_subject)) .sample_items(rng, I, per_subject) else 1:I
    for (t in items_s) {
    for (combo in .product_indices(within_sizes)) {
      level_idx <- list()
      if (length(within)) for (m in seq_along(within))
        level_idx[[within[[m]]$name]] <- combo[m]
      for (f in between) {
        n_lev <- length(f$levels)
        unit <- if (f$between == "subject") s else t
        n_unit <- if (f$between == "subject") S else I
        level_idx[[f$name]] <- ((unit - 1) * n_lev) %/% n_unit
      }
      cvals <- list(); labels <- list()
      for (f in factors) {
        li <- level_idx[[f$name]]
        labels[[f$name]] <- f$levels[li + 1]
        for (col in names(f$contrasts)) cvals[[col]] <- f$contrasts[[col]][li + 1]
      }
      rows[[length(rows) + 1L]] <- list(subject = s, item = t, labels = labels, cvals = cvals)
    }
    }
  }

  # ---- continuous predictors: one draw per unit (subject- or item-level) ----
  pred_values <- list()
  for (p in predictors) {
    unit <- p$varies_by; n_unit <- if (unit == "subject") S else I
    if (unit == "item" && !has_item) stop("predictor '", p$name, "' varies_by item but design has no items")
    pmean <- if (is.null(p$mean)) 0 else p$mean; psd <- if (is.null(p$sd)) 1 else p$sd
    vals <- numeric(n_unit)
    for (u in seq_len(n_unit)) vals[u] <- pmean + psd * rng$normal()
    pred_values[[p$name]] <- vals
  }
  if (length(predictors)) for (r_i in seq_along(rows)) for (p in predictors) {
    u <- if (p$varies_by == "subject") rows[[r_i]]$subject else rows[[r_i]]$item
    rows[[r_i]]$cvals[[p$name]] <- pred_values[[p$name]][u]
  }

  # ---- random effects (subject then item) ----
  rs <- spec$random
  b_subject <- list(); subj_cols <- character(0)
  if (!is.null(rs$subject)) {
    re <- .ranef(rs$subject); subj_cols <- re$cols
    for (s in 1:S) b_subject[[s]] <- .matvec(re$L, rng$normals(length(subj_cols)))
  }
  b_item <- list(); item_cols <- character(0)
  if (has_item && !is.null(rs$item)) {
    re <- .ranef(rs$item); item_cols <- re$cols
    for (t in 1:I) b_item[[t]] <- .matvec(re$L, rng$normals(length(item_cols)))
  }

  # ---- additional grouping factors (e.g. units nested in higher-level clusters) ----
  # Any random entry other than subject/item declares `over` (the unit it groups) and `n`
  # (the number of groups). Units are assigned to groups in equal blocks.
  extra_names <- setdiff(names(rs), c("subject", "item"))
  b_group <- list(); group_meta <- list()
  for (gname in extra_names) {
    g <- rs[[gname]]; over <- g$over; K <- g$n
    n_over <- if (over == "subject") S else I
    re <- .ranef(g)
    group_meta[[gname]] <- list(over = over, cols = re$cols,
                                group_of = ((seq_len(n_over) - 1) * K) %/% n_over)
    bg <- list(); for (gi in 0:(K - 1)) bg[[gi + 1]] <- .matvec(re$L, rng$normals(length(re$cols)))
    b_group[[gname]] <- bg
  }

  # ---- linear predictor + response (residual draws here) ----
  intercept <- spec$fixed$intercept
  coeffs <- spec$fixed$coefficients
  resp <- spec$response
  family <- resp$family; yname <- resp$name
  sigma <- resp$sigma; shift <- if (is.null(resp$shift)) 0 else resp$shift
  thresholds <- resp$thresholds; ndp <- resp$round

  n <- length(rows)
  y <- numeric(n); subj_v <- integer(n); item_v <- integer(n)
  label_cols <- vapply(factors, function(f) f$name, character(1))
  pred_names <- if (length(predictors)) vapply(predictors, function(p) p$name, character(1)) else character(0)
  label_mat <- matrix("", n, length(label_cols), dimnames = list(NULL, label_cols))
  pred_mat <- matrix(0, n, length(pred_names), dimnames = list(NULL, pred_names))
  group_mat <- matrix(0L, n, length(extra_names), dimnames = list(NULL, extra_names))

  for (r_i in seq_len(n)) {
    r <- rows[[r_i]]; cv <- r$cvals
    eta <- intercept
    for (col in names(coeffs)) eta <- eta + coeffs[[col]] * .design_value(cv, col)
    if (length(subj_cols)) {
      b <- b_subject[[r$subject]]; eta <- eta + b[1]
      if (length(subj_cols) > 1) for (j in 2:length(subj_cols))
        eta <- eta + b[j] * (if (is.null(cv[[subj_cols[j]]])) 0 else cv[[subj_cols[j]]])
    }
    if (has_item && length(item_cols)) {
      b <- b_item[[r$item]]; eta <- eta + b[1]
      if (length(item_cols) > 1) for (j in 2:length(item_cols))
        eta <- eta + b[j] * (if (is.null(cv[[item_cols[j]]])) 0 else cv[[item_cols[j]]])
    }
    for (gname in extra_names) {
      gm <- group_meta[[gname]]
      unit <- if (gm$over == "subject") r$subject else r$item
      gi <- gm$group_of[unit]; group_mat[r_i, gname] <- gi + 1L
      b <- b_group[[gname]][[gi + 1]]; eta <- eta + b[1]
      if (length(gm$cols) > 1) for (j in 2:length(gm$cols))
        eta <- eta + b[j] * (if (is.null(cv[[gm$cols[j]]])) 0 else cv[[gm$cols[j]]])
    }
    val <- switch(family,
      gaussian          = eta + sigma * rng$normal(),
      shifted_lognormal = shift + exp(eta + sigma * rng$normal()),
      lognormal         = exp(eta + sigma * rng$normal()),
      bernoulli         = if (rng$uniform() < .inv_logit(eta)) 1 else 0,
      poisson           = .poisson_inv(exp(eta), rng$uniform()),
      ordinal           = .ordinal_inv(eta, thresholds, rng$uniform()),
      beta              = { mu <- .inv_logit(eta); phi <- if (is.null(resp$phi)) 10 else resp$phi
                            .beta_draw(rng, mu * phi, (1 - mu) * phi) },
      stop("unknown family: ", family))
    if (!is.null(ndp) && family %in% c("gaussian", "shifted_lognormal", "lognormal", "beta")) val <- round(val, ndp)
    y[r_i] <- val
    subj_v[r_i] <- r$subject; item_v[r_i] <- r$item
    for (cn in label_cols) label_mat[r_i, cn] <- r$labels[[cn]]
    for (pn in pred_names) pred_mat[r_i, pn] <- cv[[pn]]
  }

  df <- data.frame(subject = subj_v, stringsAsFactors = FALSE)
  if (has_item) df$item <- item_v
  for (gname in extra_names) df[[gname]] <- group_mat[, gname]
  for (cn in label_cols) df[[cn]] <- label_mat[, cn]
  for (pn in pred_names) df[[pn]] <- pred_mat[, pn]
  df[[yname]] <- y
  df
}
