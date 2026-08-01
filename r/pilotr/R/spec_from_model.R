# Read a design specification off a model that has already been fitted.
#
# Choosing plausible random-effect standard deviations with nothing to read them from is the
# steepest step in setting up a simulation-based power analysis, and the one that most often
# stops a user before they reach the analysis itself. The tutorials on the method all say so,
# and all reach for the same remedy, which is to take the variance components from a pilot fit
# or from a published model (Green and MacLeod, 2016; Kumle et al., 2021; DeBruine and Barr,
# 2021). This file does that mechanically, mapping the estimates of a fitted lmer onto the
# fields of a pilotr specification, which can then be scaled up to the sample size being
# planned.
#
# The one thing a fitted model cannot report is which of its numeric predictor columns were
# contrast-coded factors and which were continuous covariates, because by the time the model is
# fitted both are columns of numbers. That distinction has to be inferred from the data, and it
# is inferred by a rule stated in full under spec_from_model() and reported back to the user on
# every call, because getting it wrong silently would yield a specification that simulates a
# different design while looking entirely correct.

# Refuse anything that is not a Gaussian linear mixed model, and say what to do instead.
#
# The classes worth naming separately are the ones a user is most likely to arrive with. A glmer
# fit keeps its random effects on the link scale and has no free residual standard deviation, so
# its estimates do not map onto the same fields; a brms fit and an nlme fit are not merMod
# objects at all; and a model with no random effects is missing exactly the part of a
# specification that is hardest to guess, which is the reason this function exists.
.refuse_unless_lmer <- function(fit) {
  if (inherits(fit, "lmerMod")) return(invisible(NULL))
  if (inherits(fit, "glmerMod"))
    stop("spec_from_model() reads a linear mixed model, and this is a generalised linear mixed ",
         "model fitted with glmer(). Its random effects sit on the link scale and its residual ",
         "variability is fixed by the family rather than estimated, so the estimates do not map ",
         "onto the same specification fields. Fit the model with lme4::lmer() on the scale you ",
         "mean to analyse, or write the response block of the specification by hand.",
         call. = FALSE)
  if (inherits(fit, "nlmerMod"))
    stop("spec_from_model() reads a linear mixed model, and this is a nonlinear mixed model ",
         "fitted with nlmer(), whose fixed effects are parameters of a nonlinear function ",
         "rather than coefficients on design columns.", call. = FALSE)
  hint <- if (inherits(fit, "brmsfit"))
    " For a brms fit, read the estimates off summary() and write the specification directly; brms_bridge() goes the other way, from a specification to brms code."
  else if (inherits(fit, c("lme", "nlme", "gls")))
    " For an nlme fit, refit the same model with lme4::lmer() and pass that instead."
  else if (inherits(fit, c("glm", "lm")))
    " A model with no random effects supplies no random-effect standard deviations, which are the part of a specification hardest to guess and the reason for reading one off a fit at all."
  else ""
  stop("spec_from_model() needs a linear mixed model fitted with lme4::lmer() or ",
       "lmerTest::lmer() (an object inheriting from 'lmerMod'), not an object of class ",
       paste(sprintf("'%s'", class(fit)), collapse = ", "), ".", hint, call. = FALSE)
}

# Every way `cn` could be split at its underscores into two or more parts, finest split first.
#
# model_data() writes an interaction key "a:b" as a product column "a_b", so recovering the key
# means splitting the name again. The split is ambiguous, since a column may have an underscore
# in its own name, so all the candidates are generated here and .product_parts() decides between
# them on the data. Finest first, because that is the split model_data() would have produced.
.underscore_splits <- function(cn) {
  tok <- strsplit(cn, "_", fixed = TRUE)[[1]]
  k <- length(tok)
  if (k < 2L) return(list())
  gaps <- k - 1L
  out <- vector("list", 2L^gaps - 1L)
  for (m in seq_len(2L^gaps - 1L)) {
    cut <- bitwAnd(m, bitwShiftL(1L, seq_len(gaps) - 1L)) > 0L
    parts <- character(0)
    start <- 1L
    for (i in seq_len(gaps)) if (cut[i]) {
      parts <- c(parts, paste(tok[start:i], collapse = "_"))
      start <- i + 1L
    }
    out[[m]] <- c(parts, paste(tok[start:k], collapse = "_"))
  }
  out[order(-vapply(out, length, integer(1)))]
}

# The component columns of `cn` if it really is an interaction product column, and NULL if not.
#
# The test is the product identity in the data, not the name, because a column genuinely called
# "z_freq" is a name with an underscore in it rather than an interaction between "z" and "freq".
# Verifying the arithmetic also settles the ambiguity of where to split, and it costs one pass
# over a column. The tolerance is relative, since the product of two standardised predictors is
# not reproduced to the last bit by a second multiplication.
.product_parts <- function(cn, mf, cols) {
  v <- mf[[cn]]
  if (!grepl("_", cn, fixed = TRUE) || !is.numeric(v)) return(NULL)
  scale <- max(1, max(abs(v)))
  for (parts in .underscore_splits(cn)) {
    if (!all(parts %in% cols)) next
    if (!all(vapply(parts, function(p) is.numeric(mf[[p]]), logical(1)))) next
    prod <- Reduce(`*`, lapply(parts, function(p) mf[[p]]))
    if (max(abs(prod - v)) <= 1e-8 * scale) return(parts)
  }
  NULL
}

# A fitted-model term name written in the specification's key convention.
#
# lme4 already writes an interaction between two model-frame columns as "a:b", which is what a
# specification uses, so the only work is to expand any product column back into its components
# and to rename the intercept. The expansion repeats because a product column may itself have a
# product column among its parts.
.spec_key <- function(term, products) {
  if (identical(term, "(Intercept)")) return("intercept")
  parts <- strsplit(term, ":", fixed = TRUE)[[1]]
  for (pass in seq_len(8L)) {
    if (!any(parts %in% names(products))) break
    parts <- unlist(lapply(parts, function(p)
      if (p %in% names(products)) products[[p]] else p), use.names = FALSE)
  }
  paste(parts, collapse = ":")
}

# Whether `v` takes a single value within every level of the grouping factor `g`.
.constant_within <- function(v, g) {
  all(vapply(split(v, g), function(x) length(unique(x)) <= 1L, logical(1)))
}

# The specification units within which `v` is constant, in subject-then-item order.
#
# This is what separates a unit-level variable from one that varies inside every unit, and so
# what decides `between` against `vary_within` for a factor and `varies_by` for a continuous
# predictor. Subject comes first so that a variable constant within both units, which happens in
# a nested or partially crossed design, is attributed to the subject.
.constant_units <- function(v, unit_group, flist) {
  names(unit_group)[vapply(unit_group, function(g) .constant_within(v, flist[[g]]), logical(1))]
}

# The `between` or `vary_within` field a factor should carry, given where it is constant.
#
# A factor constant within a unit assigns one level to each of that unit's members, which is
# what `between` means. A factor constant within neither unit takes both levels inside every
# unit, which is what `vary_within` means.
.placement <- function(const, unit_group) {
  if (length(const)) list(between = const[[1]]) else list(vary_within = names(unit_group))
}

# One value of `v` per level of `g`, for a unit-level variable whose moments are wanted per unit
# rather than per row. Averaging over rows would weight a unit by how many rows it contributed,
# which distorts the mean and the standard deviation of an unbalanced design.
.unit_values <- function(v, g) {
  vals <- vapply(split(v, g), function(x) if (length(x)) x[[1L]] else NA_real_, numeric(1))
  unname(vals[!is.na(vals)])
}

# Whether the levels of `inner` each sit inside a single level of `outer`.
.group_nested_in <- function(inner, outer) {
  all(vapply(split(as.integer(outer), inner),
             function(x) length(unique(x)) <= 1L, logical(1)))
}

# Whether two grouping factors cross, meaning that neither is nested in the other. Only a
# crossed second factor can play the part of pilotr's item unit; a nested one is a coarser or
# finer grouping and belongs in an extra `random` entry with its own `over` and `n`.
.groups_crossed <- function(a, b) !.group_nested_in(a, b) && !.group_nested_in(b, a)

# A name not already taken by a model-frame column, a grouping factor or an earlier factor.
.unique_name <- function(nm, used) {
  while (nm %in% used) nm <- paste0(nm, "_")
  nm
}

# A level label for one value of a two-valued numeric column.
#
# The values themselves are the only information the fit carries about what the two levels were,
# so they are also the most honest labels. .shortest_double() gives the shortest decimal that
# reads back as the same double, which keeps "-0.5" from becoming "-0.500000".
.level_label <- function(z) {
  s <- .shortest_double(z)
  if (is.null(s)) format(z) else s
}

# One grouping factor's random-effect entry, merged across the VarCorr blocks that belong to it.
#
# A `||` term splits one grouping factor across several VarCorr blocks, named "g", "g.1" and so
# on, and the split is exactly the statement that those terms were not allowed to correlate.
# Merging the blocks and setting `correlated = FALSE` records that, rather than losing it and
# letting a later model_formula() re-estimate a correlation the fit deliberately fixed at zero.
# A correlation is read only from a single undivided block, where one was actually estimated.
.random_entry <- function(vc, blocks, products) {
  sds <- numeric(0)
  cors <- list()
  for (b in blocks) {
    m <- vc[[b]]
    sd_b <- attr(m, "stddev")
    nm <- vapply(names(sd_b), .spec_key, character(1), products = products)
    sds[nm] <- as.numeric(sd_b)
    if (length(blocks) == 1L && length(nm) > 1L) {
      R <- attr(m, "correlation")
      for (i in seq_len(length(nm) - 1L)) for (j in (i + 1L):length(nm)) {
        # A correlation involving a term estimated at zero variance is NaN, and there is no
        # correlation to record in that case.
        if (is.finite(R[i, j])) cors[[paste0(nm[i], ",", nm[j])]] <- as.numeric(R[i, j])
      }
    }
  }
  # A random-effect term absent from the fit, which is what `(0 + x | g)` produces, becomes a
  # standard deviation of zero rather than a missing field, since `intercept_sd` is required and
  # zero is how a specification holds a term at no variance while keeping the structure intact.
  entry <- list(intercept_sd = if ("intercept" %in% names(sds)) sds[["intercept"]] else 0)
  slope_names <- setdiff(names(sds), "intercept")
  if (length(slope_names))
    entry$slopes <- stats::setNames(as.list(unname(sds[slope_names])), slope_names)
  if (length(cors)) entry$correlations <- cors
  if (length(blocks) > 1L && length(slope_names)) entry$correlated <- FALSE
  entry
}

# The response block, from the residual standard deviation of the fit and the family asked for.
#
# A lmer fit is Gaussian on the scale it was fitted on, so that is the default. A user who fitted
# log reaction times will want one of the lognormal families instead, and those carry parameters
# that a Gaussian fit knows nothing about, which is why `family` also accepts a list supplying
# them. Anything still missing is named in the error rather than filled in with a guess, because
# a fabricated shift or precision would change the simulated data without saying so.
.model_response <- function(fit, family, ndp, yname) {
  extra <- list()
  if (is.list(family)) {
    extra <- family[setdiff(names(family), "family")]
    family <- family$family
  }
  fam <- if (is.null(family)) "gaussian" else family
  if (!.is_scalar_string(fam) || !fam %in% names(.family_params))
    stop("`family` must name one of ", paste(names(.family_params), collapse = ", "),
         ", or be a list such as list(family = \"shifted_lognormal\", shift = 200) that also ",
         "supplies that family's parameters.", call. = FALSE)
  resp <- list(family = fam, name = yname)
  needed <- .family_params[[fam]]
  if ("sigma" %in% needed) resp$sigma <- as.numeric(stats::sigma(fit))
  for (k in names(extra)) resp[[k]] <- extra[[k]]
  absent <- setdiff(needed, names(resp))
  if (length(absent))
    stop("the ", fam, " family needs ", paste(sprintf("'%s'", absent), collapse = ", "),
         ", which cannot be read off a linear mixed model. Supply ",
         if (length(absent) > 1) "them" else "it", " with the list form of `family`, as in ",
         "family = list(family = \"", fam, "\", ",
         paste(sprintf("%s = ...", absent), collapse = ", "), ").", call. = FALSE)
  if (!is.null(ndp)) {
    if (!.is_whole(ndp) || ndp < 0)
      stop("`round` must be a single whole number of at least 0, or NULL", call. = FALSE)
    if (!fam %in% .rounding_families)
      stop("`round` has no effect for the ", fam, " family, whose outcome is already an integer",
           call. = FALSE)
    resp$round <- as.integer(ndp)
  }
  resp
}

#' Build a design specification from a fitted mixed model
#'
#' Read a design specification off a linear mixed model already fitted with `lme4::lmer()` or
#' `lmerTest::lmer()`. The fixed effects, the random-effect standard deviations and
#' correlations, and the residual standard deviation are taken from the fit, and the numbers of
#' subjects and items may be raised at the same time, so that a pilot study or a published model
#' becomes the starting point of a power analysis rather than a set of numbers to invent.
#' Requires the `lme4` package.
#'
#' @details
#' Settling on plausible random-effect standard deviations with nothing to read them from is the
#' step that most often stops a simulation-based power analysis before it starts, and the
#' tutorials on the method converge on the same remedy, which is to take the variance components
#' from a pilot fit or from a published model (Green and MacLeod, 2016; Kumle et al., 2021;
#' DeBruine and Barr, 2021). This function performs that transfer, and the specification it
#' returns can be enlarged and passed straight to [simulate_design()] or
#' [power_mixed()].
#'
#' Every number is taken on the scale the model was fitted on, and the returned family is
#' Gaussian by default, because that is what a `lmer` fit is. A model of log reaction times
#' therefore yields an intercept, coefficients and residual standard deviation on the log scale,
#' which is the right thing for the `lognormal` and `shifted_lognormal` families, whose linear
#' predictor lives on that scale as well. Ask for one of those with `family`, in the list form
#' when the family needs a parameter the fit cannot supply, as in
#' `family = list(family = "shifted_lognormal", shift = 200)`.
#'
#' Whether a numeric column of the model frame was a contrast-coded factor or a continuous
#' covariate is not recorded anywhere in a fitted model, so it is inferred here. A column that is
#' a factor, a character vector or a logical vector is treated as categorical and its own
#' contrast coding is read with [stats::contrasts()], which makes the emitted contrast-column
#' names agree with the coefficient names lme4 produced. A numeric column with exactly two
#' distinct values is also treated as categorical, on the reasoning that a two-valued numeric
#' predictor in a factorial experiment is a coded factor far more often than it is a covariate.
#' Any other numeric column becomes a continuous predictor, with `mean` and `sd` taken from the
#' data, one value per unit for a unit-level variable so that an unbalanced design does not
#' distort them. This is a heuristic and nothing more, it is reported by a message on every call,
#' and a genuine two-valued covariate has to be moved from `factors` to `predictors` by hand.
#'
#' A column is placed as `between` a unit when it holds one value within every member of that
#' unit, and as `vary_within` when it takes several values inside every unit; the same test gives
#' a continuous predictor its `varies_by`, which becomes `"observation"` when the predictor
#' varies inside every unit. A model that codes one factor twice, once as a factor for its fixed
#' effect and once as a numeric contrast for a random slope, gives two separate specification
#' terms, because the two columns are separate columns in the model frame and nothing in the fit
#' ties them together.
#'
#' Interactions need one further step. lme4 writes the interaction of two model-frame columns as
#' `a:b`, which is already the specification's convention, but [model_data()] gives an
#' interaction its own product column named `a_b`, so a specification built from a fit of
#' pilotr's own modelling data would otherwise acquire a spurious independent term. Such a column
#' is recognised by checking the product identity in the data rather than by reading its name,
#' which both avoids mistaking a column called `z_freq` for an interaction and settles where to
#' split a name with several underscores. A recognised product column is reported and re-keyed to
#' `a:b`, and it contributes no factor or predictor of its own.
#'
#' Grouping factors named `subject` and `item` are used as they stand. Otherwise the one with the
#' most levels becomes `subject`, the largest remaining factor that genuinely crosses it becomes
#' `item`, and every other grouping factor becomes an extra `random` entry with the `over` and
#' `n` fields that pilotr's additional grouping factors take, `over` being decided by which unit
#' the factor partitions. Any renaming is reported by a message and recorded in the
#' `group_mapping` attribute of the result. When each subject saw only some of the items, the
#' item unit gains a `per_subject` count, so that the recovered design keeps the partial crossing
#' of the original rather than silently becoming fully crossed.
#'
#' A boundary-singular pilot fit is carried across as it stands, which means a variance estimated
#' at zero or a correlation estimated at exactly plus or minus one. Those are faithful readings
#' of the fit rather than defects, but they are also the sign that the random-effect structure was
#' richer than the pilot could support, and a power analysis resting on them will inherit that
#' (Bates et al., 2015). Widening the design, or simplifying the structure before refitting, is
#' the remedy.
#'
#' @param fit A fitted linear mixed model, of class `lmerMod` (from `lme4::lmer()`) or
#'   `lmerModLmerTest` (from `lmerTest::lmer()`). Anything else is refused with a message saying
#'   what to pass instead.
#' @param name A label for the returned specification. Defaults to `"from_model"`.
#' @param seed The master seed of the returned specification. Defaults to 1.
#' @param n_subject Number of subjects for the returned specification. Defaults to the number of
#'   levels the fit actually had, and is normally raised above it, since scaling a pilot design
#'   up is the point of reading a specification off a pilot fit.
#' @param n_item Number of items, treated the same way as `n_subject`. An error when the model
#'   has no second crossed grouping factor to act as items, because there is then no item unit to
#'   resize.
#' @param family The response family of the returned specification. `NULL`, the default, gives
#'   `"gaussian"`, which is what a `lmer` fit is on the scale it was fitted on. A single string
#'   names another family, and a list such as
#'   `list(family = "shifted_lognormal", shift = 200)` also supplies the parameters of that
#'   family which a Gaussian fit cannot provide.
#' @param round Decimal places for the simulated response, passed through to `response.round`.
#'   `NULL`, the default, leaves the response unrounded.
#' @return A design specification as a nested list, carrying `spec_version` and validated with
#'   [validate_spec()], so it can be passed directly to [simulate_design()],
#'   [power_mixed()] or [spec_json()]. Two attributes record the readings that
#'   the fit did not settle on its own: `group_mapping`, a named character vector giving the
#'   specification unit each of the fit's grouping factors became, and `column_kinds`, a named
#'   character vector giving each model-frame column the classification `"factor"` or
#'   `"predictor"`. Both are attributes rather than fields so that the specification itself stays
#'   within the portable schema.
#' @references Bates, D., Kliegl, R., Vasishth, S. and Baayen, H. (2015). Parsimonious mixed
#'   models. \emph{arXiv}. \doi{10.48550/arXiv.1506.04967}
#'
#'   DeBruine, L. M. and Barr, D. J. (2021). Understanding mixed-effects models through data
#'   simulation. \emph{Advances in Methods and Practices in Psychological Science}, 4(1).
#'   \doi{10.1177/2515245920965119}
#'
#'   Green, P. and MacLeod, C. J. (2016). SIMR: an R package for power analysis of generalized
#'   linear mixed models by simulation. \emph{Methods in Ecology and Evolution}, 7(4), 493-498.
#'   \doi{10.1111/2041-210X.12504}
#'
#'   Kumle, L., Vo, M. L.-H. and Draschkow, D. (2021). Estimating power in (generalized) linear
#'   mixed models: An open introduction and tutorial in R. \emph{Behavior Research Methods}, 53,
#'   2528-2543. \doi{10.3758/s13428-021-01546-0}
#' @seealso [simulate_design()] to simulate from the recovered specification,
#'   [power_mixed()] to run the power analysis it was read off the fit for, and
#'   [model_formula()] for the analysis model a specification implies.
#' @examples
#' \donttest{
#' if (requireNamespace("lme4", quietly = TRUE)) {
#'   # Stand in for a pilot study: simulate a small design and fit it.
#'   pilot <- build_spec(list(name = "pilot", seed = 1, design_kind = "within",
#'     include_items = TRUE, n_subject = 30, n_item = 20, factor_name = "cond",
#'     lev1 = "a", lev2 = "b", intercept = 6, effect = 0.05,
#'     subj_int_sd = 0.12, subj_slope_sd = 0, subj_corr = 0,
#'     item_int_sd = 0.08, item_slope_sd = 0, item_corr = 0,
#'     family = "shifted_lognormal", resp_name = "", sigma = 0.3, shift = 200))
#'   fit <- lme4::lmer(model_formula(pilot), data = model_data(pilot, simulate_design(pilot)))
#'
#'   # Read the design back off the fit, scaled up to the sample size being planned.
#'   spec <- spec_from_model(fit, n_subject = 60, n_item = 40)
#'   spec$random$subject$intercept_sd
#'   attr(spec, "column_kinds")
#'   head(simulate_design(spec))
#' }
#' }
#' @export
spec_from_model <- function(fit, name = NULL, seed = 1, n_subject = NULL, n_item = NULL,
                            family = NULL, round = NULL) {
  if (!requireNamespace("lme4", quietly = TRUE))
    stop("spec_from_model() requires the 'lme4' package; please install it.", call. = FALSE)
  .refuse_unless_lmer(fit)
  if (!is.null(name) && !.is_scalar_string(name))
    stop("`name` must be a single string, or NULL", call. = FALSE)
  if (!.is_whole(seed)) stop("`seed` must be a single whole number", call. = FALSE)
  for (a in list(list("n_subject", n_subject), list("n_item", n_item)))
    if (!is.null(a[[2]]) && (!.is_whole(a[[2]]) || a[[2]] < 1))
      stop("`", a[[1]], "` must be a single whole number of at least 1, or NULL", call. = FALSE)

  mf <- stats::model.frame(fit)
  flist <- lme4::getME(fit, "flist")
  gnames <- unique(names(flist))

  # ---- the response column ----
  # The formula's left-hand side names the column the specification should generate, so that the
  # recovered specification feeds back into model_data() and model_formula() unchanged. A
  # transformed outcome such as log(rt) is not a column name, so its variable is used instead.
  lhs <- stats::formula(fit)[[2]]
  yname <- paste(deparse(lhs), collapse = "")
  if (!grepl("^[.]?[A-Za-z][A-Za-z0-9._]*$", yname)) {
    vars <- all.vars(lhs)
    yname <- if (length(vars)) vars[[1L]] else "outcome"
  }

  # ---- grouping factors mapped onto the specification's units ----
  nlev <- vapply(gnames, function(g) nlevels(flist[[g]]), integer(1))
  by_size <- gnames[order(-nlev)]
  subj <- if ("subject" %in% gnames) "subject" else by_size[[1L]]
  rest <- setdiff(by_size, subj)
  itm <- NULL
  if ("item" %in% rest) itm <- "item"
  else for (g in rest) if (.groups_crossed(flist[[subj]], flist[[g]])) { itm <- g; break }
  extras <- setdiff(rest, itm)
  fit_groups <- c(subj, itm, extras)
  mapping <- stats::setNames(c("subject", if (!is.null(itm)) "item", extras), fit_groups)
  renamed <- names(mapping)[names(mapping) != mapping]
  if (length(renamed))
    message("spec_from_model() mapped the grouping factor",
            if (length(renamed) > 1L) "s " else " ",
            paste(sprintf("'%s' onto '%s'", renamed, mapping[renamed]), collapse = ", "),
            ", because a specification names its units 'subject' and 'item'.")
  unit_group <- c(subject = subj)
  if (!is.null(itm)) unit_group["item"] <- itm

  # ---- units ----
  units <- list(subject = list(
    n = as.integer(if (is.null(n_subject)) nlev[[subj]] else n_subject)))
  if (!is.null(itm)) {
    units$item <- list(n = as.integer(if (is.null(n_item)) nlev[[itm]] else n_item))
    # Partial crossing, where each subject saw a subset of the items, changes the number of rows
    # and the precision that follows from it, so it is worth carrying over rather than quietly
    # returning a fully crossed design. The count is clamped to the item total in case the caller
    # asked for fewer items than the pilot had.
    seen <- vapply(split(as.integer(flist[[itm]]), flist[[subj]]),
                   function(x) length(unique(x)), integer(1))
    if (length(seen) && max(seen) < nlev[[itm]]) {
      per <- max(1L, min(as.integer(stats::median(seen)), units$item$n))
      units$item$per_subject <- per
      message("spec_from_model() found each subject seeing about ", per, " of the ", nlev[[itm]],
              " items, so the item unit carries per_subject = ", per, ".")
    }
  } else if (!is.null(n_item)) {
    stop("`n_item` was supplied, but this model has no second crossed grouping factor to treat ",
         "as items, so there is no item unit to resize.", call. = FALSE)
  }

  # ---- design columns, and the product columns among them ----
  drop <- unique(c(yname, all.vars(lhs), gnames))
  cols <- setdiff(names(mf), drop)
  cols <- cols[!grepl("^[(]", cols)]   # (weights) and (offset) are not design columns
  products <- list()
  for (cn in cols) {
    parts <- .product_parts(cn, mf, cols)
    if (!is.null(parts)) products[[cn]] <- parts
  }
  design_cols <- setdiff(cols, names(products))

  # ---- factors and continuous predictors ----
  factors <- list()
  predictors <- list()
  kinds <- character(0)
  used <- c(names(mf), gnames)
  for (cn in design_cols) {
    v <- mf[[cn]]
    const <- .constant_units(v, unit_group, flist)
    if (is.factor(v) || is.character(v) || is.logical(v)) {
      f <- if (is.factor(v)) v else factor(v)
      cm <- stats::contrasts(f)
      lab <- colnames(cm)
      if (is.null(lab) || !all(nzchar(lab))) lab <- as.character(seq_len(ncol(cm)))
      contr <- stats::setNames(
        lapply(seq_len(ncol(cm)), function(j) as.numeric(cm[, j])), paste0(cn, lab))
      factors[[length(factors) + 1L]] <- c(
        list(name = cn, levels = levels(f), contrasts = contr), .placement(const, unit_group))
      used <- c(used, names(contr))
      kinds[cn] <- "factor"
    } else if (!is.numeric(v)) {
      stop("column '", cn, "' of the model frame is of class '",
           paste(class(v), collapse = "/"), "', which a specification has no field for. ",
           "Recode it as a factor or a numeric column and refit.", call. = FALSE)
    } else if (length(unique(v)) == 2L) {
      vals <- sort(unique(v))
      fname <- .unique_name(paste0(cn, "_level"), used)
      used <- c(used, fname)
      factors[[length(factors) + 1L]] <- c(
        list(name = fname, levels = vapply(vals, .level_label, character(1)),
             contrasts = stats::setNames(list(as.numeric(vals)), cn)),
        .placement(const, unit_group))
      kinds[cn] <- "factor"
    } else {
      vb <- if (length(const)) const[[1L]] else "observation"
      vals <- if (identical(vb, "observation")) v else .unit_values(v, flist[[unit_group[[vb]]]])
      predictors[[length(predictors) + 1L]] <- list(
        name = cn, varies_by = vb, mean = mean(vals),
        sd = if (length(vals) > 1L) stats::sd(vals) else 0)
      kinds[cn] <- "predictor"
    }
  }
  facs <- names(kinds)[kinds == "factor"]
  preds <- names(kinds)[kinds == "predictor"]
  if (length(kinds))
    message("spec_from_model() read ",
            if (length(facs)) paste0(paste(sprintf("'%s'", facs), collapse = ", "),
                                     " as categorical ",
                                     if (length(facs) > 1L) "factors" else "factor") else "",
            if (length(facs) && length(preds)) " and " else "",
            if (length(preds)) paste0(paste(sprintf("'%s'", preds), collapse = ", "),
                                      " as continuous ",
                                      if (length(preds) > 1L) "predictors" else "predictor") else "",
            ". A fitted model does not record which of its numeric columns were coded factors, ",
            "so please check that reading.")
  if (length(products))
    message("spec_from_model() read ",
            paste(sprintf("'%s' as the interaction '%s'", names(products),
                          vapply(names(products), .spec_key, character(1),
                                 products = products)), collapse = ", "), ".")

  # ---- fixed effects ----
  fe <- lme4::fixef(fit)
  keys <- setdiff(names(fe), "(Intercept)")
  coefs <- stats::setNames(lapply(keys, function(k) as.numeric(fe[[k]])),
                           vapply(keys, .spec_key, character(1), products = products))
  # A model fitted without an intercept still needs the field, and zero is what its absence
  # means. Saying so is better than leaving a required field out and failing validation.
  if (!"(Intercept)" %in% names(fe))
    message("spec_from_model() found no intercept in the fit, so fixed.intercept is 0.")
  icept <- if ("(Intercept)" %in% names(fe)) as.numeric(fe[["(Intercept)"]]) else 0

  # ---- random effects ----
  vc <- lme4::VarCorr(fit)
  base_of <- function(nm) if (nm %in% gnames) nm else sub("[.][0-9]+$", "", nm)
  vc_base <- vapply(names(vc), base_of, character(1))
  random <- list()
  for (fg in fit_groups) {
    entry <- .random_entry(vc, names(vc)[vc_base == fg], products)
    sg <- mapping[[fg]]
    if (!sg %in% c("subject", "item")) {
      entry$over <- if (.group_nested_in(flist[[subj]], flist[[fg]])) "subject"
        else if (!is.null(itm) && .group_nested_in(flist[[itm]], flist[[fg]])) "item"
        else {
          message("spec_from_model() could not tell which unit '", fg, "' groups, so it is ",
                  "recorded as grouping subjects. Check random$", sg, "$over.")
          "subject"
        }
      entry$n <- as.integer(nlev[[fg]])
    }
    random[[sg]] <- entry
  }

  # ---- assemble ----
  spec <- list(spec_version = .SPEC_VERSION,
               name = if (is.null(name)) "from_model" else name,
               seed = as.integer(seed),
               units = units)
  if (length(factors)) spec$factors <- factors
  if (length(predictors)) spec$predictors <- predictors
  spec$fixed <- list(intercept = icept, coefficients = coefs)
  spec$random <- random
  spec$response <- .model_response(fit, family, round, yname)
  validate_spec(spec)
  structure(spec, group_mapping = mapping, column_kinds = kinds)
}
