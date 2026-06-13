# Spec parsing + generative engine -- a bit-identical mirror of simdgp/simulate.py.

#' Load a JSON design specification.
#' @export
load_spec <- function(path) {
  jsonlite::fromJSON(path, simplifyVector = TRUE, simplifyDataFrame = FALSE,
                     simplifyMatrix = FALSE)
}

# itertools.product order: last index varies fastest. Returns a list of 0-based index vectors.
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

# (columns, lower-Cholesky L) for one unit's random-effect covariance.
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

#' Simulate a data set from a design specification (path or parsed list).
#' @export
simulate_design <- function(spec) {
  if (is.character(spec)) spec <- load_spec(spec)

  S <- spec$units$subject$n
  has_item <- !is.null(spec$units$item)
  I <- if (has_item) spec$units$item$n else 1L

  factors <- spec$factors
  within <- Filter(function(f) !is.null(f$vary_within), factors)
  between <- Filter(function(f) !is.null(f$between), factors)
  within_sizes <- vapply(within, function(f) length(f$levels), integer(1))

  # ---- canonical row order (no randomness) ----
  rows <- list()
  for (s in 1:S) for (t in 1:I) {
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

  # ---- random effects in documented order ----
  rng <- make_rng(spec$seed)
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

  # ---- linear predictor + response (residual draws here) ----
  intercept <- spec$fixed$intercept
  coeffs <- spec$fixed$coefficients
  resp <- spec$response
  family <- resp$family; yname <- resp$name
  sigma <- resp$sigma; shift <- if (is.null(resp$shift)) 0 else resp$shift
  thresholds <- resp$thresholds; ndp <- resp$round

  n <- length(rows)
  y <- numeric(n)
  subj_v <- integer(n); item_v <- integer(n)
  label_cols <- vapply(factors, function(f) f$name, character(1))
  label_mat <- matrix("", n, length(label_cols), dimnames = list(NULL, label_cols))

  for (r_i in seq_len(n)) {
    r <- rows[[r_i]]
    eta <- intercept
    for (col in names(coeffs)) eta <- eta + coeffs[[col]] * (if (is.null(r$cvals[[col]])) 0 else r$cvals[[col]])
    if (length(subj_cols)) {
      b <- b_subject[[r$subject]]; eta <- eta + b[1]
      if (length(subj_cols) > 1) for (j in 2:length(subj_cols))
        eta <- eta + b[j] * (if (is.null(r$cvals[[subj_cols[j]]])) 0 else r$cvals[[subj_cols[j]]])
    }
    if (has_item && length(item_cols)) {
      b <- b_item[[r$item]]; eta <- eta + b[1]
      if (length(item_cols) > 1) for (j in 2:length(item_cols))
        eta <- eta + b[j] * (if (is.null(r$cvals[[item_cols[j]]])) 0 else r$cvals[[item_cols[j]]])
    }
    val <- switch(family,
      gaussian          = eta + sigma * rng$normal(),
      shifted_lognormal = shift + exp(eta + sigma * rng$normal()),
      bernoulli         = if (rng$uniform() < .inv_logit(eta)) 1 else 0,
      poisson           = .poisson_inv(exp(eta), rng$uniform()),
      ordinal           = .ordinal_inv(eta, thresholds, rng$uniform()),
      stop("unknown family: ", family))
    if (!is.null(ndp) && family %in% c("gaussian", "shifted_lognormal")) val <- round(val, ndp)
    y[r_i] <- val
    subj_v[r_i] <- r$subject; item_v[r_i] <- r$item
    for (cn in label_cols) label_mat[r_i, cn] <- r$labels[[cn]]
  }

  df <- data.frame(subject = subj_v, stringsAsFactors = FALSE)
  if (has_item) df$item <- item_v
  for (cn in label_cols) df[[cn]] <- label_mat[, cn]
  df[[yname]] <- y
  df
}
