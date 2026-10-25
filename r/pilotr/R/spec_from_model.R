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
# The fixed effects are mapped too. A pilot's point estimates are a poor basis for planning,
# since they are noisy and, when a pilot is followed up because it looked promising, biased
# upwards (Albers and Lakens, 2018). Every call therefore says which effects it took from the
# fit, and the caller can replace them with safeguard values or a smallest effect size of interest.
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

# The component columns of `cn` if it is an interaction product column, and NULL if not.
#
# The test is the product identity in the data, because a column called "z_freq" may be nothing
# more than a name with an underscore in it, with no "z" and "freq" behind it.
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

# Whether a factor placed within both units was counterbalanced in the pilot, meaning that every
# observed subject-item pair holds one level and every subject and every item holds two or more.
#
# A priming study shows each subject each item once and rotates the levels across lists of
# subjects, the counterbalanced design of Judd, Westfall and Kenny (2017). A `vary_within` factor
# crosses every subject-item pair with every level, so the specification read off such a pilot
# observes each pair at every level, which multiplies the observations and overstates power. A
# crossed pilot keeps both levels in most pairs even when rows are missing, so one pair with two
# levels is enough to rule the rotation out. Requiring two levels in each subject and each item
# keeps a pilot in which some subjects or items met a single level from being described as a
# rotation across lists.
.counterbalanced <- function(v, subj, item) {
  n_levels <- function(g) vapply(split(v, g, drop = TRUE),
                                 function(x) length(unique(x)), integer(1))
  pair <- paste(as.integer(subj), as.integer(item))
  all(n_levels(pair) == 1L) && all(n_levels(subj) >= 2L) && all(n_levels(item) >= 2L)
}

# One value of `v` per level of `g`, for a unit-level variable whose moments are wanted per
# unit. Averaging over rows would weight a unit by how many rows it contributed,
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

# The factors placed between the same unit, merged into a single factor of cells.
#
# Up to specification version 0.3, pilotr assigned the levels of each between factor to blocks of
# its unit on its own, so two factors read as between subjects fell into the same or overlapping
# blocks. A balanced 2 x 2 between-subjects pilot came back as cells of 40, 0, 0 and 40 at 80
# subjects, and its refit dropped the second effect and the interaction. Version 0.4 allocates
# such factors jointly. A factor of cells keeps every cell in the same way, in a specification
# that the releases reading 0.3 read too. Its contrast columns are the components' columns under
# their fitted names, so every coefficient key the fit produced still names a column, an
# interaction such as Aa2:Bb2 included. The cells run in .product_indices() order, the
# first-listed factor varying slowest, and their labels join the components' labels with ".".
.combine_between <- function(factors, used) {
  for (unit in c("subject", "item")) {
    at <- which(vapply(factors, function(f) identical(f[["between"]], unit), logical(1)))
    if (length(at) < 2L) next
    parts <- factors[at]
    nms <- vapply(parts, function(f) f[["name"]], character(1))
    cells <- .product_indices(vapply(parts, function(f) length(f[["levels"]]), integer(1)))
    labels <- vapply(cells, function(ix) paste(vapply(seq_along(parts), function(k)
      parts[[k]][["levels"]][ix[k] + 1L], character(1)), collapse = "."), character(1))
    contr <- list()
    for (k in seq_along(parts)) for (cn in names(parts[[k]][["contrasts"]]))
      contr[[cn]] <- vapply(cells, function(ix) parts[[k]][["contrasts"]][[cn]][ix[k] + 1L],
                            numeric(1))
    name <- .unique_name(paste(nms, collapse = "_"), used)
    used <- c(used, name)
    message("spec_from_model() combined the factors ", .name_list(sprintf("'%s'", nms)), ", ",
            if (length(nms) == 2L) "both" else "all", " constant within ", unit,
            ", into one factor '", name, "' whose levels are their cells. It allocates the ",
            unit, "s as the separate factors would under specification version 0.4, and ",
            "releases that read only version 0.3 read it too.")
    # make.unique() matters only when a label itself contains ".", which could give two cells
    # one label.
    factors[[at[1L]]] <- list(name = name, levels = make.unique(labels), contrasts = contr,
                              between = unit)
    factors <- factors[-at[-1L]]
  }
  factors
}

# One grouping factor's random-effect entry, merged across the VarCorr blocks that belong to it.
#
# A `||` term splits one grouping factor across several VarCorr blocks, named "g", "g.1" and so
# on, and the split is exactly the statement that those terms were not allowed to correlate.
# Merging the blocks and setting `correlated = FALSE` records that. Losing it would let a later
# model_formula() re-estimate a correlation the fit deliberately fixed at zero.
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
  # standard deviation of zero, since `intercept_sd` is required and
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
# them. Anything still missing is named in the error, because a guessed shift or precision would
# change the simulated data without saying so.
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

# Three significant digits, enough for a message to identify a value. The returned specification
# holds the full number.
.sig3 <- function(x) vapply(x, function(v) format(v, digits = 3), character(1), USE.NAMES = FALSE)

# The `sesoi` argument that a refusal suggests for the named effects, each value left as a
# placeholder. A key that is not a syntactic name, such as the interaction "a:b", is quoted with
# backticks, as c() needs it.
.sesoi_example <- function(keys) {
  quoted <- ifelse(make.names(keys) == keys, keys, paste0("`", keys, "`"))
  paste0("sesoi = c(", paste(sprintf("%s = ...", quoted), collapse = ", "), ")")
}

# Refuse a `sesoi` that is not a set of named numbers. Whether the names are coefficients of the
# returned specification is settled once the fit's terms have been given the specification's keys.
.check_sesoi <- function(sesoi) {
  if (is.null(sesoi)) return(invisible(NULL))
  nm <- names(sesoi)
  if (!is.numeric(sesoi) || !length(sesoi) || is.null(nm) || anyNA(nm) || !all(nzchar(nm)) ||
      anyDuplicated(nm) > 0L || !all(is.finite(sesoi)))
    stop("`sesoi` must be NULL or a numeric vector of finite values with distinct, non-empty ",
         "names, as in sesoi = c(cond = 0.03)", call. = FALSE)
  invisible(NULL)
}

# The bound nearer zero of each named term's two-sided interval at `level`, which is the safeguard
# value of Perugini, Gallucci and Costantini (2014), with the interval itself and whether its
# quantiles came from a t distribution.
#
# The interval is the estimate plus or minus a quantile times the standard error from vcov(fit).
# For an lme4 fit, the quantile is normal, which treats the standard error as known, so the
# interval is too narrow when the subjects or items are few. A lmerTest fit carries
# Satterthwaite's degrees of freedom for each term, and a t quantile on them allows for the
# uncertainty in the variance components, although only approximately. Unless lmerTest supplies
# them for every term, every term takes the normal quantile, so that one message describes them
# all.
.safeguard_bounds <- function(fit, terms, level) {
  est <- lme4::fixef(fit)[terms]
  se <- sqrt(diag(as.matrix(stats::vcov(fit))))[terms]
  df <- NULL
  if (inherits(fit, "lmerModLmerTest") && requireNamespace("lmerTest", quietly = TRUE)) {
    co <- tryCatch(suppressMessages(summary(fit, ddf = "Satterthwaite"))$coefficients,
                   error = function(e) NULL)
    if (!is.null(co) && "df" %in% colnames(co) && all(terms %in% rownames(co)) &&
        all(is.finite(co[terms, "df"]) & co[terms, "df"] > 0))
      df <- co[terms, "df"]
  }
  p <- (1 + level) / 2
  q <- if (is.null(df)) stats::qnorm(p) else stats::qt(p, df)
  lo <- est - q * se
  hi <- est + q * se
  list(lo = lo, hi = hi, value = ifelse(est > 0, lo, hi), t = !is.null(df))
}

# The message saying where each coefficient of the returned specification came from, whether
# `sesoi`, the fit's estimate or a safeguard value read off the fit. When any is the pilot's own
# estimate, the message also says why that is a poor value to plan a study for.
.effects_message <- function(coefs, estimates, set, method, level, t_dist) {
  listed <- function(k) .name_list(sprintf("'%s' = %s", k, .sig3(unlist(coefs[k]))))
  taken <- setdiff(names(coefs), set)
  one <- length(taken) == 1L
  clauses <- if (length(set)) paste0("set ", listed(set), " from `sesoi`") else character(0)
  if (length(taken))
    clauses <- c(clauses, if (method == "estimate")
      paste0("took ", listed(taken), " from the fit, as the pilot estimated ",
             if (one) "it" else "them")
    else
      paste0("took ", listed(taken), " from the fit as ",
             if (one) "a safeguard value" else "safeguard values", ", the bound nearer zero of ",
             if (one) "its " else "each one's ", format(100 * level), "% interval, computed with ",
             if (t_dist) "the t distribution with Satterthwaite's degrees of freedom"
             else "the normal distribution"))
  msg <- paste0("spec_from_model() ", paste(clauses, collapse = " and "), ".")
  if (!length(taken)) return(msg)
  if (method == "safeguard")
    return(paste0(msg, " The pilot estimated ", if (one) "it" else "them", " at ",
                  .name_list(.sig3(estimates[taken])), "."))
  paste0(msg, " A pilot's estimates are noisy, and a study powered for them tends to be ",
         "underpowered. Set a smallest effect size of interest with `sesoi`, or use ",
         "effects = \"safeguard\" (see ?spec_from_model).")
}

#' Build a design specification from a fitted mixed model
#'
#' Read a design specification off a linear mixed model already fitted with `lme4::lmer()` or
#' `lmerTest::lmer()`. The fixed effects, the random-effect standard deviations and
#' correlations, and the residual standard deviation are taken from the fit, and the numbers of
#' subjects and items may be raised at the same time, so that a pilot study or a published model
#' becomes the starting point of a power analysis. The fixed effects of a pilot are its point
#' estimates, and a study powered for them tends to be underpowered. A message names them on
#' every call, and `sesoi` or `effects = "safeguard"` replaces them (see Details).
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
#' The fixed effects are another matter. They are the pilot's point estimates, which a small
#' pilot measures imprecisely, and a pilot followed up because its effect looked promising has
#' on average overestimated that effect. A study powered for such an estimate therefore tends to
#' have less power than planned (Albers and Lakens, 2018). The variance components are the part
#' of a pilot worth carrying over. The effect is better set to a smallest effect size of interest
#' (SESOI), the smallest effect of theoretical or practical importance. For an effect keyed
#' `cond`, as in the `crossed_mixed_rt` example, pass `sesoi = c(cond = 0.03)`, or set
#' `spec$fixed$coefficients$cond <- 0.03` on the returned specification. Every call reports in a
#' message which coefficients it took from the fit, with their values.
#'
#' When no SESOI can be argued for, `effects = "safeguard"` takes safeguard values (Perugini,
#' Gallucci and Costantini, 2014). Each coefficient not named in `sesoi` becomes the bound of its
#' two-sided `safeguard_level` interval nearer zero, by default that of the 60% interval. The
#' interval is the estimate plus or minus a quantile times the standard error from `vcov(fit)`.
#' For an `lme4` fit, the quantile is normal, `qnorm((1 + safeguard_level) / 2)`, which is 0.84 at
#' the default. The normal quantile treats the standard error as known, so the interval is too
#' narrow when the subjects or items are few. For a `lmerTest` fit, the quantile is a t quantile
#' on the Satterthwaite degrees of freedom that `lmerTest` reports, which allow for the
#' uncertainty in the variance components, although only approximately. When `lmerTest` cannot
#' give them for every coefficient, all coefficients take the normal quantile, and the message
#' names the distribution used. An interval that includes zero stops the call with an error
#' naming the effect. The pilot then leaves even the sign of the effect open, and a value of zero
#' would plan a study of no effect at all. Such an effect needs a value in `sesoi`. By default,
#' the conditions that [design_conditions()] builds for [sweep_spec()] include one with every
#' named effect at zero, which shows how often a design declares an effect that is not there.
#'
#' Whatever their source, the coefficients of the specification are the effects its data are
#' simulated with. A named numeric `focal` in [power_mixed()] changes only the true values against
#' which Type S and Type M are computed. The simulated effect stays as it was, so powering a study
#' for another effect means changing the specification.
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
#' Two or more factors placed between the same unit become one factor whose levels are their
#' cells, and a message says so. The first-listed factor varies slowest, and each label joins the
#' components' labels with `"."`. Two between factors over one unit need specification version
#' 0.4, which allocates them jointly. Up to version 0.3, each was assigned to blocks of its unit on
#' its own, which confounded their effects, and [validate_spec()] refuses them in a specification
#' declaring an earlier version. The factor of cells allocates the units as version 0.4 allocates
#' the separate factors, and every release that reads version 0.3 reads it, so the returned
#' specification declares 0.3. The combined factor carries every component's contrast columns
#' under the names the fit gave them, so each coefficient the fit estimated, interactions
#' included, still names a column. Every cell then receives units, in equal numbers whenever the
#' number of units is a multiple of the number of cells.
#'
#' Interactions need one further step. lme4 writes the interaction of two model-frame columns as
#' `a:b`, which is already the specification's convention, but [model_data()] gives an
#' interaction its own product column named `a_b`, so a specification built from a fit of
#' pilotr's own modelling data would otherwise acquire a spurious independent term. Such a column
#' is recognised by checking the product identity in the data, which both avoids mistaking a
#' column called `z_freq` for an interaction and settles where to
#' split a name with several underscores. A recognised product column is reported and re-keyed to
#' `a:b`, and it contributes no factor or predictor of its own.
#'
#' Grouping factors named `subject` and `item` are used as they stand. Otherwise the one with the
#' most levels becomes `subject`, the largest remaining factor that is crossed with it becomes
#' `item`, and every other grouping factor becomes an extra `random` entry with the `over` and
#' `n` fields that pilotr's additional grouping factors take, `over` being decided by which unit
#' the factor partitions. Any renaming is reported by a message and recorded in the
#' `group_mapping` attribute of the result. When each subject saw only some of the items, the
#' item unit gains a `per_subject` count, so that the recovered design keeps the partial crossing
#' of the original.
#'
#' A factor that varies within both units is crossed with every subject-item pair unless it sets
#' `counterbalance`, a field of specification version 0.4. A pilot run in counterbalanced lists,
#' as priming studies usually are, showed each subject each item at one level only. The levels
#' rotated across lists of subjects, so that every item appeared at each level. When every
#' subject-item pair of the pilot holds one level of a factor, and every subject and every item
#' holds two or more, the returned specification sets `counterbalance` to `TRUE` on that factor
#' and declares version 0.4. Each subject then meets each item at one level, the levels rotating
#' over subjects and items as described under Counterbalanced lists in the specification
#' (`spec/SPEC.md`). The simulated data name each subject's list in a column `list`. Factors
#' rotated together, as two two-level factors over four lists, are all counterbalanced, and a
#' message names them. Releases that read no version above 0.3 refuse such a specification.
#'
#' `counterbalance = FALSE` returns the crossed specification instead, with a warning, and releases
#' that read no version above 0.3 accept it. It has twice the pilot's observations per subject for a
#' two-level factor, three times for a three-level one and four times for two two-level factors
#' rotated together. The power computed from it is therefore too high. Take the variance
#' components of the `crossed_mixed_rt` example, with 30 subjects, 24 items and an effect of 0.05
#' on the log scale. By the formulas of Judd, Westfall and Kenny (2017), power is then about 0.80
#' for the crossed design and 0.55 for the counterbalanced one. The counterbalanced study needs
#' about 1.8 times the subjects to match the crossed one. A release that reads no version above 0.3
#' takes a two-level design written with a list factor between subjects and an item-set factor
#' between items, as in the two-list example among the worked encodings of the specification.
#' The crossed specification and the warning also come back, whatever `counterbalance` says, when
#' each subject of the pilot saw only some of the items. The specification then keeps
#' `per_subject`, which version 0.4 does not combine with `counterbalance`, since that field needs
#' every subject to see every item. They come back as well when the pilot's model has a column or
#' grouping factor named `list`, as a model that fits the list beside the condition does. Version
#' 0.4 gives that name to the column holding each subject's list, and no other column may share
#' it. The warning then asks for the term to be renamed in the pilot's data and the model refitted.
#'
#' Alternatively, keep the crossed specification and drop the rows that a counterbalanced study
#' would not observe, through the `prep` argument of [power_mixed()] and [precision_design()].
#' That function receives each simulated data set and returns the modelling data, so it ends by
#' calling [model_data()]. With `workers > 1` it runs in worker processes that have not attached
#' pilotr and cannot see the caller's variables, so it has to call `pilotr::model_data()` and
#' carry the specification with it. A function made by another function carries it once that
#' function has forced its argument, as in this rotation over two lists for the
#' `crossed_mixed_rt` example, whose factor is `condition`:
#'
#' ```r
#' counterbalance <- function(spec) {
#'   force(spec)
#'   function(d) {
#'     level <- match(d$condition, spec$factors[[1]]$levels) - 1L
#'     pilotr::model_data(spec, d[(d$subject %% 2L + d$item + level) %% 2L == 0L, ])
#'   }
#' }
#' power_mixed(spec, prep = counterbalance(spec), workers = 4)
#' ```
#'
#' A boundary-singular pilot fit is carried across as it stands, which means a variance estimated
#' at zero or a correlation estimated at exactly plus or minus one. Those are faithful readings
#' of the fit, and no kind of defect, but they are also the sign that the random-effect
#' structure was richer than the pilot could support, and a power analysis resting on them will
#' inherit that
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
#' @param effects Where the coefficients of the returned specification come from. `"estimate"`,
#'   the default, takes the fit's point estimates. `"safeguard"` takes
#'   the bound of each coefficient's two-sided `safeguard_level` interval nearer zero, and stops
#'   with an error when that interval includes zero (see Details). Either way, a coefficient named
#'   in `sesoi` takes the value given there, and the intercept is the fit's.
#' @param safeguard_level The confidence level of the intervals behind `effects = "safeguard"`, a
#'   single number greater than 0 and less than 1. The default of 0.6 follows Perugini, Gallucci
#'   and Costantini (2014).
#' @param sesoi `NULL`, the default, or a named numeric vector that sets coefficients to a smallest
#'   effect size of interest on the scale of the fit, as in `sesoi = c(cond = 0.03)`. The names
#'   are the returned specification's coefficient keys, which write an interaction as `a:b`, and
#'   a name that is not among them is refused.
#' @param counterbalance Whether a factor that the pilot counterbalanced sets `counterbalance` in
#'   the returned specification, which then declares specification version 0.4 (see Details).
#'   `TRUE`, the default, plans the pilot's design, unless its subjects each saw only some of the
#'   items or its model has a term named `list`. `FALSE` returns a fully crossed specification
#'   that releases reading version 0.3 accept, with a warning that the power computed from it is
#'   too high.
#' @return A design specification as a nested list, carrying the lowest `spec_version` that reads
#'   it, 0.4 when it sets `counterbalance` and 0.3 otherwise, and validated with
#'   [validate_spec()], so it can be passed directly to [simulate_design()],
#'   [power_mixed()] or [spec_json()]. Two attributes record the readings that
#'   the fit did not settle on its own: `group_mapping`, a named character vector giving the
#'   specification unit each of the fit's grouping factors became, and `column_kinds`, a named
#'   character vector giving each model-frame column the classification `"factor"` or
#'   `"predictor"`. A third, `effects`, is a list recording where the coefficients came from. Its
#'   `method` is the value of `effects`, its `level` the `safeguard_level` applied, or `NA` with
#'   `"estimate"`, and its `sesoi` the keys of the coefficients set from `sesoi`. All three are
#'   attributes rather than fields so that the specification itself stays within the portable
#'   schema.
#' @references Albers, C. and Lakens, D. (2018). When power analyses based on pilot data are
#'   biased: Inaccurate effect size estimators and follow-up bias. \emph{Journal of Experimental
#'   Social Psychology}, 74, 187-195. \doi{10.1016/j.jesp.2017.09.004}
#'
#'   Bates, D., Kliegl, R., Vasishth, S. and Baayen, H. (2015). Parsimonious mixed
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
#'   Judd, C. M., Westfall, J. and Kenny, D. A. (2017). Experiments with more than one random
#'   factor: Designs, analytic models, and statistical power. \emph{Annual Review of Psychology},
#'   68, 601-625. \doi{10.1146/annurev-psych-122414-033702}
#'
#'   Kumle, L., Vo, M. L.-H. and Draschkow, D. (2021). Estimating power in (generalized) linear
#'   mixed models: An open introduction and tutorial in R. \emph{Behavior Research Methods}, 53,
#'   2528-2543. \doi{10.3758/s13428-021-01546-0}
#'
#'   Perugini, M., Gallucci, M. and Costantini, G. (2014). Safeguard power as a protection against
#'   imprecise power estimates. \emph{Perspectives on Psychological Science}, 9(3), 319-332.
#'   \doi{10.1177/1745691614528519}
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
#'
#'   # The effect is the pilot's estimate. Plan for a smallest effect size of interest, here 0.03
#'   # on the log scale, or for the bound of the estimate's 60% interval nearer zero.
#'   spec_from_model(fit, n_subject = 60, n_item = 40, sesoi = c(effect = 0.03))$fixed
#'   spec_from_model(fit, n_subject = 60, n_item = 40, effects = "safeguard")$fixed
#' }
#' }
#' @export
spec_from_model <- function(fit, name = NULL, seed = 1, n_subject = NULL, n_item = NULL,
                            family = NULL, round = NULL, effects = c("estimate", "safeguard"),
                            safeguard_level = 0.6, sesoi = NULL, counterbalance = TRUE) {
  if (!requireNamespace("lme4", quietly = TRUE))
    stop("spec_from_model() requires the 'lme4' package; please install it.", call. = FALSE)
  .refuse_unless_lmer(fit)
  if (!is.null(name) && !.is_scalar_string(name))
    stop("`name` must be a single string, or NULL", call. = FALSE)
  if (!.is_whole(seed)) stop("`seed` must be a single whole number", call. = FALSE)
  for (a in list(list("n_subject", n_subject), list("n_item", n_item)))
    if (!is.null(a[[2]]) && (!.is_whole(a[[2]]) || a[[2]] < 1))
      stop("`", a[[1]], "` must be a single whole number of at least 1, or NULL", call. = FALSE)
  effects <- match.arg(effects)
  if (!.is_scalar_number(safeguard_level) || safeguard_level <= 0 || safeguard_level >= 1)
    stop("`safeguard_level` must be a single number greater than 0 and less than 1",
         call. = FALSE)
  .check_sesoi(sesoi)
  if (!isTRUE(counterbalance) && !isFALSE(counterbalance))
    stop("`counterbalance` must be TRUE or FALSE", call. = FALSE)

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
  # The columns of factors that the pilot counterbalanced, with each factor's number of levels.
  # A crossed specification multiplies the observations per subject by their product.
  # `rotated_factor` gives the factor each column became, which takes another name when the
  # column holds two numbers, as a -0.5/0.5 column `cond` becomes the factor `cond_level`.
  rotated <- integer(0)
  rotated_factor <- character(0)
  is_rotated <- function(v, const)
    !length(const) && !is.null(itm) && .counterbalanced(v, flist[[subj]], flist[[itm]])
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
      if (is_rotated(v, const)) {
        rotated[cn] <- nlevels(f)
        rotated_factor[cn] <- cn
      }
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
      if (is_rotated(v, const)) {
        rotated[cn] <- 2L
        rotated_factor[cn] <- fname
      }
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
  factors <- .combine_between(factors, c(used, yname))

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

  # ---- where the coefficients come from ----
  # They are the pilot's point estimates unless the caller chose otherwise. `sesoi` is keyed as the
  # returned specification is, so an interaction is "a:b" even when the fit named it after the
  # product column "a_b".
  spec_keys <- as.character(names(coefs))
  unknown <- setdiff(names(sesoi), spec_keys)
  if (length(unknown))
    stop("`sesoi` names ", .name_list(sprintf("'%s'", unknown)),
         if (length(unknown) == 1L) ", which is not a coefficient"
         else ", which are not coefficients", " of the returned specification.",
         if (length(spec_keys)) paste0(" Its coefficients are ",
                                       .name_list(sprintf("'%s'", spec_keys)), ".")
         else " It has none.", call. = FALSE)
  set <- spec_keys[spec_keys %in% names(sesoi)]
  for (k in set) coefs[[k]] <- as.numeric(sesoi[[k]])
  estimates <- stats::setNames(vapply(keys, function(k) as.numeric(fe[[k]]), numeric(1)),
                               spec_keys)
  rest <- setdiff(spec_keys, set)
  t_dist <- FALSE
  if (effects == "safeguard" && length(rest)) {
    b <- .safeguard_bounds(fit, keys[match(rest, spec_keys)], safeguard_level)
    # An interval that holds zero leaves even the sign of the effect open. Its bound nearer zero
    # is zero itself, which would plan a study of no effect at all, so the caller is asked for a
    # value instead.
    open <- !(is.finite(b$lo) & is.finite(b$hi) & (b$lo > 0 | b$hi < 0))
    if (any(open)) {
      ko <- rest[open]
      one <- length(ko) == 1L
      ivals <- sprintf("from %s to %s", .sig3(b$lo[open]), .sig3(b$hi[open]))
      stop("spec_from_model() cannot give ", .name_list(sprintf("'%s'", ko)),
           if (one) " a safeguard value: its " else " safeguard values: their ",
           format(100 * safeguard_level), "% interval",
           if (one) paste0(", ", ivals, ", includes zero")
           else paste0("s include zero (", .name_list(sprintf("%s for '%s'", ivals, ko)), ")"),
           ", so the pilot does not settle ",
           if (one) "the sign of the effect" else "the signs of these effects", ". Give ",
           if (one) "it a smallest effect size" else "them smallest effect sizes",
           " of interest with `sesoi`, as in ", .sesoi_example(ko), ".", call. = FALSE)
    }
    for (i in seq_along(rest)) coefs[[rest[i]]] <- unname(b$value[[i]])
    t_dist <- b$t
  }
  if (length(coefs))
    message(.effects_message(coefs, estimates, set, effects, safeguard_level, t_dist))

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

  # ---- counterbalanced factors ----
  # Specification version 0.4 writes the pilot's design with `counterbalance` on each factor the
  # pilot rotated, so that the specification shows each subject each item once, as the pilot did.
  # The field needs every subject to see every item, so a pilot whose subjects each saw some of the
  # items, which keeps per_subject, stays crossed, as does every pilot when the caller asks for a
  # specification that the releases reading 0.3 accept. The field also gives the data a column
  # `list`, which no other column may share. A pilot whose model has a term of that name, as when
  # the list is fitted beside the condition, stays crossed too, since validation would refuse it.
  partial <- !is.null(units[["item"]][["per_subject"]])
  list_taken <- "list" %in% c(names(random), yname,
                              vapply(factors, function(f) f[["name"]], character(1)),
                              unlist(lapply(factors, function(f) names(f[["contrasts"]]))),
                              vapply(predictors, function(p) p[["name"]], character(1)))
  set_field <- length(rotated) > 0L && counterbalance && !partial && !list_taken
  if (set_field)
    for (k in seq_along(factors))
      if (factors[[k]][["name"]] %in% rotated_factor) factors[[k]][["counterbalance"]] <- TRUE

  # ---- assemble ----
  spec <- list(name = if (is.null(name)) "from_model" else name,
               seed = .as_seed(seed),
               units = units)
  if (length(factors)) spec$factors <- factors
  if (length(predictors)) spec$predictors <- predictors
  spec$fixed <- list(intercept = icept, coefficients = coefs)
  spec$random <- random
  spec$response <- .model_response(fit, family, round, yname)
  # The lowest version that reads the specification, first, as spec_json() writes it. The factor
  # of cells built above keeps a factorial between design within 0.3, where two between factors
  # would need 0.4, so the releases that read 0.3 read this specification too, unless it sets
  # counterbalance, a field of 0.4.
  spec <- c(list(spec_version = .required_spec_version(spec)), spec)
  validate_spec(spec)
  # The pilot names its columns, and the specification its factors, which differ for a column of
  # two numbers. Factors rotated together, as over the four lists of a 2 x 2 design, are named
  # once, with the product of their level counts.
  if (length(rotated)) {
    one <- length(rotated) == 1L
    times <- prod(rotated)
    every <- if (one) "level" else "combination"
    on <- .name_list(sprintf("'%s'", rotated_factor))
    seen <- paste0(
      "spec_from_model(): in the pilot each subject saw each item under one ",
      if (one) sprintf("level of '%s'", names(rotated))
      else paste("combination of the levels of", .name_list(sprintf("'%s'", names(rotated)))),
      " (a counterbalanced, Latin-square design)")
    if (set_field) {
      message(seen, ". The returned specification therefore sets counterbalance to true on ", on,
              " and declares spec_version \"0.4\", so that each subject meets each item at one ",
              every, ", as in the pilot (see Counterbalanced lists in SPEC.md). Releases that ",
              "read no version above 0.3 refuse it. With counterbalance = FALSE, ",
              "spec_from_model() returns a crossed specification that they read, but the power ",
              "computed from it is too high.")
    } else {
      # A warning, where the readings above are messages, because this specification would plan
      # a larger study than the pilot ran. The power computed from it would be too high, and
      # nothing downstream would show it. For one factor of two levels, the warning also names the
      # two-list encoding, which releases reading 0.3 take and which allows per_subject. Where the
      # pilot's model takes the name `list`, setting the field by hand would be refused, so the
      # advice is to rename that term, whichever `counterbalance` was asked for.
      two_list <- paste0("list (between subject) x item set (between item), as in SPEC.md's ",
                         "two-list example.")
      warning(seen, ", but the returned specification ",
              if (partial) paste0("crosses the items each subject sees with every ", every)
              else paste0("crosses every subject with every item under every ", every), ", ",
              if (times == 2) "doubling the observations per subject"
              else paste("multiplying the observations per subject by", times),
              " and overstating power.",
              if (partial)
                paste0(" Specification version 0.4 counterbalances a factor only when every ",
                       "subject sees every item, and in the pilot each subject saw some of the ",
                       "items, so the returned specification keeps units.item.per_subject and ",
                       "does not set counterbalance.",
                       if (times == 2) paste0(" Encode the design as ", two_list))
              else if (list_taken)
                paste0(" Specification version 0.4 names each subject's list in a column 'list', ",
                       "a name that the pilot's model uses already, so ",
                       if (counterbalance) "the returned specification does not set counterbalance"
                       else "a specification setting counterbalance would be refused",
                       ". To plan the pilot's design, rename 'list' in the pilot's data, refit ",
                       "the model and call spec_from_model() again",
                       if (!counterbalance) " with counterbalance = TRUE", ".",
                       if (times == 2) paste0(" For a release that reads no version above 0.3, ",
                                              "encode it as ", two_list))
              else
                paste0(" To plan the pilot's design, set counterbalance to true on ", on,
                       " and declare spec_version \"0.4\" (see Counterbalanced lists in SPEC.md).",
                       if (times == 2) paste0(" For a release that reads no version above 0.3, ",
                                              "encode it as ", two_list)),
              call. = FALSE)
    }
  }
  structure(spec, group_mapping = mapping, column_kinds = kinds,
            effects = list(method = effects,
                           level = if (effects == "safeguard") safeguard_level else NA_real_,
                           sesoi = set))
}
