# A sample-size justification, written in Markdown from a power result.
#
# The documentation calls the value target_n() returns the number a preregistration quotes, yet
# the justification that goes with it had to be assembled by hand. It names the analysis and its
# test, each effect size and where it came from, the solved sample size with its interval and the
# power when the null hypothesis is true. A simulation adds its Monte Carlo error and fit counts,
# and a reader needs what reproduces the run. power_report() writes all of these down.
#
# It works in two steps. .power_report_record() reads a result, a solved curve and a specification
# file into a record of counts and strings, and .report_render() writes the record as Markdown.
# Every rate is computed in the renderer from two counts, the significant replicates and the
# replicates behind them, with .rate_with_error() and .type_s_interval() from mcse.R. The Python
# twin, python/pilotr/report.py, renders a record with the same operations and the same format
# strings, so one record gives the same bytes in both languages, and tools/parity/report_cross.py
# checks this over the records in report_cases.json. Numbers go through sprintf() here and the %
# operator there. Both round a double correctly to the digits asked for, an exact tie going to the
# even digit, so they write the same text for the same double.
#
# The section on intervals and the region of practical equivalence, for a precision_design()
# result, exists in R alone, as precision_design() does.

#' Write a sample-size justification from a power result
#'
#' Turns the result of [power_mixed()], [power_design()] or [precision_design()], with the
#' [target_n()] result that solved a curve for the sample size, into a sample-size justification in
#' Markdown, ready for a preregistration, a registered report or a grant application.
#'
#' @details
#' The sections follow the recommendations of Lakens (2022, Table 4) for reporting an a-priori
#' power analysis. The primary analysis comes first, with its model, test, alpha and sidedness,
#' then the desired power and the sample size solved for it, with its interval. Each focal effect
#' follows, with its true value, the justification given for that value and its power. Last comes
#' the power when the null hypothesis is true. The simulation is then accounted for as Morris, White
#' and Crowther (2019) recommend for simulation studies: the replicates simulated and those that
#' returned an estimate, the fits that were singular or carried a warning, and the Monte Carlo
#' standard error of every rate. The report ends with what a reader needs to reproduce the run:
#' the software and its versions, the specification's seed and an MD5 checksum of the
#' specification file's bytes. The versions are those loaded when the report is written, so write
#' it in the session that ran the analysis.
#'
#' Each rate's Monte Carlo standard error and 95% Wilson interval are computed in the report from
#' the number of significant replicates and the number of replicates behind the rate, with the
#' functions that compute them in the results. Type S has its Wilson interval from the number of
#' wrong signs among the significant replicates. The Python twin's `power_report()` writes the
#' same text, byte for byte, from the same counts, which `tools/parity/report_cross.py` in the
#' source repository checks.
#'
#' The rates are those of the design `x` simulated. A sample size solved from a curve need not be
#' that design's, so a justification usually pairs the solved size with a run of [power_mixed()]
#' at that size. `spec_file` names the specification of that run, and the report lists its units.
#' A result carries neither its specification nor its seed, so without `spec_file` the report says
#' that both are missing.
#'
#' `null_result` is the result of the same analysis of a design whose focal effects are all zero,
#' such as the null condition that [design_conditions()] puts first. Its rejection rate is the
#' Type I error rate of the test, which for a Wald z with few subjects or items can exceed alpha,
#' and which a preregistered simplification rule can raise (see [power_mixed()]).
#'
#' A [precision_design()] result gives a report on the interval and the region of practical
#' equivalence in place of the test. Each focal effect then has its mean interval width, the
#' probability of a meaningful effect, with the interval wholly outside the region, and the
#' probability of practical equivalence, with the interval wholly inside it. `rope` states the
#' region, since the result does not record it.
#' A [target_n()] result solved over `p_meaningful` gives the desired probability of a meaningful
#' decision. This section exists in R alone, as [precision_design()] does.
#'
#' @param x For `power_report()`, the result of [power_mixed()], [power_design()] or
#'   [precision_design()] at one design. A curve goes to [target_n()], and its result to `solved`.
#'   For `print()`, a report.
#' @param solved Optional result of [target_n()], or of [solve_curve()] for an axis other than
#'   sample size, which gives the desired power and the value solved for it.
#' @param spec_file Optional path of the specification file the analysis read. The report gives
#'   the MD5 checksum of its bytes, its seed and its units.
#' @param justification Optional justification of each focal effect's size: one string for every
#'   focal effect, or a character vector named by focal effect. An effect without one is reported
#'   as having none.
#' @param null_result Optional result of the same analysis with every focal effect set to zero.
#'   A result whose focal effects are not zero is refused.
#' @param rope The region of practical equivalence given to [precision_design()], as a half-width
#'   or an interval, which the report states. It applies to a precision result alone.
#' @param ... Ignored, present for consistency with the generic.
#' @return A character string of class `pilotr_report` holding the Markdown, which ends in a
#'   newline. Printing it shows the text, and `cat(report, file = "justification.md")` saves it.
#'   The same arguments always give the same text.
#' @references Lakens, D. (2022). Sample size justification. \emph{Collabra: Psychology}, 8(1),
#'   33267. \doi{10.1525/collabra.33267}
#'
#'   Morris, T. P., White, I. R. and Crowther, M. J. (2019). Using simulation studies to evaluate
#'   statistical methods. \emph{Statistics in Medicine}, 38(11), 2074-2102.
#'   \doi{10.1002/sim.8086}
#' @examples
#' spec <- build_spec(list(name = "d", seed = 1, design_kind = "between",
#'   factor_name = "group", lev1 = "a", lev2 = "b", n_subject = 64,
#'   intercept = 100, effect = 5, family = "gaussian", resp_name = "", sigma = 10))
#' # n_sims is small so the example runs quickly. Use 200 or more for real planning.
#' res <- power_design(spec, n_sims = 50)
#' spec_file <- tempfile(fileext = ".json")
#' writeLines(spec_json(spec), spec_file)
#' power_report(res, spec_file = spec_file,
#'   justification = "Half a standard deviation, the smallest difference worth acting on.")
#' @seealso [target_n()] for the sample size a curve implies, and [design_conditions()] for the
#'   null condition.
#' @export
power_report <- function(x, solved = NULL, spec_file = NULL, justification = NULL,
                         null_result = NULL, rope = NULL) {
  rec <- .power_report_record(x, solved = solved, spec_file = spec_file,
                              justification = justification, null_result = null_result,
                              rope = rope)
  structure(.report_render(rec), class = c("pilotr_report", "character"))
}

#' @rdname power_report
#' @export
print.pilotr_report <- function(x, ...) {
  cat(x)
  invisible(x)
}

# The refusals. Those about `solved`, `justification`, `null_result` and `spec_file` are word for
# word those of the Python twin.
.REPORT_X <- paste0(
  "`x` must be the result of power_mixed(), power_design() or precision_design() at one design. ",
  "A curve goes to target_n(), and its result to `solved`.")
.REPORT_SOLVED <- "`solved` must be the result of target_n() or solve_curve()."
.REPORT_JUSTIFICATION <- paste0(
  "`justification` must be one string, which applies to every focal effect, or strings named by ",
  "focal effect.")
.REPORT_UNNAMED <- "`x` has one effect, without a name, so `justification` must be a single string."
.REPORT_NULL_KIND <- "`null_result` must be a result of the same kind as `x`."
.REPORT_ROPE <- "`rope` states the region of practical equivalence of a precision_design() result."
.REPORT_SPEC_FILE <- "`spec_file` must be the path of the specification file the analysis read."

# The heading of the one effect of a power_design() result, which names none.
.REPORT_TWO_GROUP <- "The difference between the two groups"

# What a result is: "power" for power_mixed() or power_design() and "precision" for
# precision_design(), each at one design. A curve, from a sweep, is refused, since its rows are
# designs of their own.
.report_kind <- function(x) {
  if (is.data.frame(x)) {
    if (length(x) && names(x)[1] == "param" && "p_meaningful" %in% names(x) && nrow(x))
      return("precision")
    stop(.REPORT_X, call. = FALSE)
  }
  if (is.list(x) && all(c("power", "n_significant", "n_sims", "alpha") %in% names(x)) &&
      length(x$power) >= 1L && (inherits(x, "pilotr_power") || length(x$power) == 1L))
    return("power")
  stop(.REPORT_X, call. = FALSE)
}

# Element i of a per-effect vector, or NULL when it is missing.
.report_value <- function(v, i) {
  if (is.null(v) || length(v) < i || is.na(v[[i]])) NULL else unname(v[[i]])
}

# The number of replicates behind each effect's power. power_mixed() takes an effect's power over
# the fits that gave it an estimate, which can be fewer than n_returned when a coefficient is
# absent from some fits, and the result keeps the rate and the count of significant replicates.
# The denominator is their ratio, which division recovers exactly. With no significant replicate
# it is n_returned, and power_design() takes power over every replicate.
.report_tested <- function(p, k, base) {
  if (is.na(p)) 0 else if (p > 0) round(k / p) else base
}

.report_power_effects <- function(x) {
  nm <- names(x$power)
  base <- if (!is.null(x$n_returned)) x$n_returned else x$n_sims
  lapply(seq_along(x$power), function(i) {
    name <- if (is.null(nm)) NULL else nm[[i]]
    k <- unname(x$n_significant[[i]])
    list(name = name,
         label = if (is.null(name)) .REPORT_TWO_GROUP else NULL,
         coefficients = if (!is.null(x$sets) && !is.null(name)) as.list(x$sets[[name]]) else NULL,
         true_effect = .report_value(x$true_effect, i),
         n_significant = k,
         n_tested = .report_tested(unname(x$power[[i]]), k, base),
         type_s = .report_value(x$type_s, i),
         type_m = .report_value(x$type_m, i),
         type_m_mcse = .report_value(x$type_m_mcse, i))
  })
}

# A precision_design() result keeps each decision as a proportion over the replicates that
# returned an estimate, and the count is that proportion times their number.
.report_precision_effects <- function(x) {
  lapply(seq_len(nrow(x)), function(i) {
    n <- x$n_returned[[i]]
    count <- function(p) if (is.na(p)) NULL else round(p * n)
    list(name = x$param[[i]], label = NULL, coefficients = NULL,
         true_effect = .report_value(x$true, i),
         n_tested = n,
         n_meaningful = count(x$p_meaningful[[i]]),
         n_equivalent = count(x$p_equivalent[[i]]),
         mean_ci_width = .report_value(x$mean_ci_width, i))
  })
}

# Each focal effect's justification, or NULL where none was given.
.report_justify <- function(justification, effects) {
  if (is.null(justification)) return(rep(list(NULL), length(effects)))
  if (!is.character(justification) || !length(justification) || anyNA(justification))
    stop(.REPORT_JUSTIFICATION, call. = FALSE)
  given <- names(justification)
  if (is.null(given)) {
    if (length(justification) != 1L) stop(.REPORT_JUSTIFICATION, call. = FALSE)
    return(rep(list(unname(justification[[1]])), length(effects)))
  }
  if (any(!nzchar(given)) || anyDuplicated(given)) stop(.REPORT_JUSTIFICATION, call. = FALSE)
  names_ <- vapply(effects, function(e) if (is.null(e$name)) NA_character_ else e$name,
                   character(1))
  if (anyNA(names_)) stop(.REPORT_UNNAMED, call. = FALSE)
  unknown <- setdiff(given, names_)
  if (length(unknown))
    stop(sprintf("`justification` names '%s', which is not a focal effect of `x`. Its focal effects are %s.",
                 unknown[1], .name_list(sprintf("'%s'", names_))), call. = FALSE)
  lapply(names_, function(n) if (n %in% given) unname(justification[[n]]) else NULL)
}

# The fitter a power result names. power_design() names none and always fits a pooled-variance t.
.report_fitter <- function(x) {
  if (!is.null(x$fitter)) return(x$fitter)
  if (!inherits(x, "pilotr_power")) return("stats::t.test, two-sample t with pooled variance")
  NULL
}

# power_design() records its alternative. power_mixed() names a one-sided test in its fitter, and
# earlier results were two-sided.
.report_alternative <- function(x) {
  if (!is.null(x$alternative)) return(x$alternative)
  f <- x$fitter
  if (is.null(f) || !grepl("one-sided", f, fixed = TRUE)) return("two.sided")
  if (grepl("greater", f, fixed = TRUE)) "greater" else "less"
}

# The share of the replicates analysed at each rung of a ladder, from power_mixed()'s rung_shares
# or from precision_design()'s rung_<name> columns.
.report_rungs <- function(x, kind) {
  shares <- if (kind == "power") x$rung_shares else {
    cols <- grep("^rung_", names(x), value = TRUE)
    if (length(cols)) stats::setNames(vapply(cols, function(k) x[[k]][[1]], numeric(1)),
                                      sub("^rung_", "", cols))
  }
  if (is.null(shares)) return(NULL)
  lapply(seq_along(shares), function(i) list(name = names(shares)[i], share = unname(shares[[i]])))
}

.report_fits <- function(x, kind) {
  first <- function(k) if (is.null(x[[k]])) NULL else x[[k]][[1]]
  n_sims <- if (kind == "power") x$n_sims else first("n_attempted")
  list(n_sims = n_sims, n_attempted = first("n_attempted"), n_returned = first("n_returned"),
       n_converged = first("n_converged"), n_singular = first("n_singular"),
       n_warning = first("n_warning"))
}

# The null condition's counts, after checking that its focal effects are zero.
.report_null <- function(null_result, kind) {
  if (is.null(null_result)) return(NULL)
  nk <- tryCatch(.report_kind(null_result), error = function(e) NA_character_)
  if (!identical(nk, kind)) stop(.REPORT_NULL_KIND, call. = FALSE)
  effects <- if (kind == "power") .report_power_effects(null_result)
             else .report_precision_effects(null_result)
  for (e in effects) {
    v <- e$true_effect
    if (!is.null(v) && v != 0)
      stop(sprintf("`null_result` must come from a design in which every focal effect is zero. Its effect '%s' has a true value of %s.",
                   if (is.null(e$name)) e$label else e$name, .report_g(v)), call. = FALSE)
  }
  fits <- .report_fits(null_result, kind)
  keep <- if (kind == "power") c("name", "label", "n_significant", "n_tested")
          else c("name", "label", "n_meaningful", "n_tested")
  list(alpha = if (kind == "power") null_result$alpha else NULL,
       alternative = if (kind == "power") .report_alternative(null_result) else NULL,
       n_sims = fits$n_sims, n_returned = fits$n_returned,
       effects = lapply(effects, function(e) e[keep]))
}

.report_solved <- function(solved) {
  if (is.null(solved)) return(NULL)
  need <- c("value", "lo", "hi", "target", "level", "transform", "x", "y", "dispersion")
  if (!is.list(solved) || !all(need %in% names(solved))) stop(.REPORT_SOLVED, call. = FALSE)
  keep <- c("target", "level", "value", "lo", "hi", "n", "n_lo", "n_hi", "transform",
            "dispersion", "x", "y", "n_points", "x_min", "x_max")
  stats::setNames(lapply(keep, function(k) solved[[k]]), keep)
}

# "pilotr 0.3.1, R 4.6.1", with lme4 and lmerTest when the fitter used them.
.report_software <- function(fitter) {
  version <- function(pkg)
    if (requireNamespace(pkg, quietly = TRUE)) paste(pkg, utils::packageVersion(pkg)) else ""
  f <- if (is.null(fitter)) "" else fitter
  parts <- c(paste("pilotr", utils::packageVersion("pilotr")), paste("R", getRversion()))
  if (grepl("lme4", f, fixed = TRUE) || grepl("lmerTest", f, fixed = TRUE))
    parts <- c(parts, version("lme4"))
  if (grepl("lmerTest", f, fixed = TRUE)) parts <- c(parts, version("lmerTest"))
  paste(parts[nzchar(parts)], collapse = ", ")
}

# The specification file's name, the MD5 checksum of its bytes, its seed and its units. A result
# records none of them. The checksum is taken over the file's bytes as they stand, with no
# decoding, as the Python twin's hashlib takes it, so the same file gives the same checksum in both.
.report_reproducibility <- function(spec_file, fitter) {
  out <- list(software = .report_software(fitter), spec_file = NULL, spec_md5 = NULL,
              seed = NULL, units = NULL)
  if (is.null(spec_file)) return(out)
  if (!is.character(spec_file) || length(spec_file) != 1L || is.na(spec_file))
    stop(.REPORT_SPEC_FILE, call. = FALSE)
  if (!file.exists(spec_file) || dir.exists(spec_file))
    stop(sprintf("%s No file is found at '%s'.", .REPORT_SPEC_FILE, spec_file), call. = FALSE)
  spec <- load_spec(spec_file, validate = FALSE)
  units <- spec[["units"]]
  out["spec_file"] <- list(basename(spec_file))
  out["spec_md5"] <- list(unname(tools::md5sum(spec_file)))
  out["seed"] <- list(spec[["seed"]])
  out["units"] <- list(lapply(names(units), function(u) list(name = u, n = units[[u]][["n"]])))
  out
}

# The record a report is rendered from: counts and strings, laid out as tools/parity/
# report_cases.json writes them, with NULL for anything missing.
.power_report_record <- function(x, solved, spec_file, justification, null_result, rope) {
  kind <- .report_kind(x)
  if (!is.null(rope)) {
    if (kind != "precision") stop(.REPORT_ROPE, call. = FALSE)
    rope <- .as_rope(rope)
  }
  effects <- if (kind == "power") .report_power_effects(x) else .report_precision_effects(x)
  just <- .report_justify(justification, effects)
  for (i in seq_along(effects)) effects[[i]]["justification"] <- list(just[[i]])
  fitter <- if (kind == "power") .report_fitter(x) else if (!is.null(x$fitter)) x$fitter[[1]]
  analysis <- list(fitter = fitter,
                   alpha = if (kind == "power") x$alpha else NULL,
                   alternative = if (kind == "power") .report_alternative(x) else NULL,
                   joint = kind == "power" && !is.null(x$sets),
                   rungs = .report_rungs(x, kind),
                   fallback = if (kind == "power") x$fallback else NULL,
                   rope = rope)
  list(kind = kind, analysis = analysis, effects = effects, solved = .report_solved(solved),
       null = .report_null(null_result, kind), fits = .report_fits(x, kind),
       reproducibility = .report_reproducibility(spec_file, fitter))
}


# --- Rendering. -----------------------------------------------------------------------------------
#
# Every line below has a twin in python/pilotr/report.py, except those that render a
# precision_design() result, which exist in R alone, as the top of this file says.
#
# Fields are read with [[ ]], which matches a name exactly. A record read from JSON leaves out a
# field it has no value for, and `$` would then match a longer name: `solved$n` found `n_points`
# in the result of solve_curve(), which has no `n`, and reported a sample size it did not solve.

# A missing value: NULL, as JSON's null reads, or NA, as a result holds it.
.report_missing <- function(v) {
  is.null(v) || !length(v) || (is.atomic(v) && length(v) == 1L && is.na(v))
}
.report_f3 <- function(v) sprintf("%.3f", as.numeric(v))
.report_g <- function(v) sprintf("%g", as.numeric(v))
.report_int <- function(v) sprintf("%.0f", as.numeric(v))
.report_name <- function(e) {
  if (.report_missing(e[["name"]])) e[["label"]] else paste0("`", e[["name"]], "`")
}

.REPORT_NO_ESTIMATE <- "not available, since no replicate returned an estimate for this effect"

# A rate with its Monte Carlo standard error and Wilson interval, computed from its counts.
.report_rate <- function(count, n) {
  if (.report_missing(n) || as.numeric(n) <= 0 || .report_missing(count))
    return(.REPORT_NO_ESTIMATE)
  r <- .rate_with_error(as.numeric(count), as.numeric(n), "p")
  sprintf("%s (Monte Carlo standard error %s; 95%% Wilson interval %s to %s), from %s of %s replicates",
          .report_f3(r[["p"]]), .report_f3(r[["p_mcse"]]), .report_f3(r[["p_lo"]]),
          .report_f3(r[["p_hi"]]), .report_int(count), .report_int(n))
}

.report_sidedness <- function(alternative, joint) {
  if (isTRUE(joint)) return("a joint test of each set")
  if (.report_missing(alternative) || alternative == "two.sided") return("two-sided")
  sprintf("one-sided (%s)", alternative)
}

.report_units <- function(units) {
  paste(vapply(units, function(u) paste(u[["name"]], .report_int(u[["n"]])), character(1)),
        collapse = ", ")
}

.report_section_analysis <- function(rec) {
  a <- rec[["analysis"]]
  fitter <- if (.report_missing(a[["fitter"]])) "not recorded" else a[["fitter"]]
  rope <- a[["rope"]]
  out <- c("## Primary analysis", "")
  if (identical(rec[["kind"]], "precision")) {
    region <- if (.report_missing(rope)) " given to precision_design()"
              else sprintf(", from %s to %s,", .report_g(rope[[1]]), .report_g(rope[[2]]))
    out <- c(out, paste0("- Model and interval: ", fitter),
             sprintf("- Decision rule: a 95%% Wald interval wholly outside the region of practical equivalence%s is a meaningful effect, and one wholly inside it is practical equivalence",
                     region))
  } else {
    level <- if (.report_missing(a[["alpha"]])) "not recorded"
             else sprintf("alpha = %s, %s", .report_g(a[["alpha"]]),
                          .report_sidedness(a[["alternative"]], a[["joint"]]))
    out <- c(out, paste0("- Model and test: ", fitter), paste0("- Significance level: ", level))
  }
  rungs <- a[["rungs"]]
  if (!.report_missing(rungs)) {
    rule <- paste(vapply(rungs, function(r) r[["name"]], character(1)), collapse = ", then ")
    fallback <- a[["fallback"]]
    if (!.report_missing(fallback))
      rule <- paste0(rule, ", moving on from ",
                     if (fallback == "warning") "singular fits and fits with warnings"
                     else "singular fits")
    out <- c(out, paste0("- Simplification rule: ", rule))
  }
  units <- rec[["reproducibility"]][["units"]]
  if (!.report_missing(units)) out <- c(out, paste0("- Simulated units: ", .report_units(units)))
  out
}

.report_section_solved <- function(rec) {
  s <- rec[["solved"]]
  if (.report_missing(s)) return(NULL)
  x <- s[["x"]]
  whole <- !.report_missing(s[["n"]])
  rate <- if (identical(s[["y"]], "p_meaningful")) "Desired probability of a meaningful decision"
          else "Desired power"
  level <- sprintf("%s%%", .report_g(as.numeric(s[["level"]]) * 100))
  out <- c(if (whole) "## Desired power and sample size" else "## Desired power and the solved value",
           "", sprintf("- %s: %s", rate, .report_g(s[["target"]])))
  out <- c(out, if (whole)
    sprintf("- Sample size: %s (%s confidence interval %s to %s), the exact solution %s rounded up",
            .report_int(s[["n"]]), level, .report_int(s[["n_lo"]]), .report_int(s[["n_hi"]]),
            sprintf("%.2f", as.numeric(s[["value"]])))
  else
    sprintf("- Solved value of `%s`: %s (%s confidence interval %s to %s)", x,
            .report_g(s[["value"]]), level, .report_g(s[["lo"]]), .report_g(s[["hi"]])))
  scale <- switch(s[["transform"]], sqrt = sprintf("the square root of `%s`", x),
                  log = sprintf("the logarithm of `%s`", x), sprintf("`%s`", x))
  span <- if (.report_missing(s[["n_points"]])) "" else
    sprintf(", over %s points from %s to %s", .report_int(s[["n_points"]]),
            .report_g(s[["x_min"]]), .report_g(s[["x_max"]]))
  c(out, sprintf("- Fit: a probit regression of `%s` on %s%s", s[["y"]], scale, span),
    sprintf("- Heterogeneity factor: %s, the dispersion that scales the variance behind the interval, which is 1 where the curve fits the probit model",
            sprintf("%.2f", as.numeric(s[["dispersion"]]))))
}

.report_justification <- function(e) {
  paste0("- Justification: ",
         if (.report_missing(e[["justification"]])) "not given" else e[["justification"]])
}

.report_effect_power <- function(e, joint) {
  out <- c(paste0("### ", .report_name(e)), "")
  if (!.report_missing(e[["coefficients"]]))
    out <- c(out, paste0("- Coefficients tested jointly: ",
                         paste(unlist(e[["coefficients"]]), collapse = ", ")))
  truth <- e[["true_effect"]]
  known <- !.report_missing(truth)
  if (!isTRUE(joint))
    out <- c(out, paste0("- True value: ", if (known) .report_g(truth)
                         else "not given, so Type S and Type M are not defined"))
  out <- c(out, .report_justification(e),
           paste0("- Power: ", .report_rate(e[["n_significant"]], e[["n_tested"]])))
  if (isTRUE(joint) || !known) return(out)
  s <- e[["type_s"]]
  type_s <- if (.report_missing(s)) "not available" else {
    k <- as.numeric(e[["n_significant"]])
    ci <- .type_s_interval(round(as.numeric(s) * k), k, as.numeric(truth))
    if (is.na(ci[1])) .report_f3(s)
    else sprintf("%s (95%% Wilson interval %s to %s)", .report_f3(s), .report_f3(ci[1]),
                 .report_f3(ci[2]))
  }
  m <- e[["type_m"]]
  m_mcse <- e[["type_m_mcse"]]
  type_m <- if (.report_missing(m)) "not available"
            else if (.report_missing(m_mcse))
              sprintf("%s (too few significant replicates for a Monte Carlo standard error)",
                      .report_f3(m))
            else sprintf("%s (Monte Carlo standard error %s)", .report_f3(m), .report_f3(m_mcse))
  c(out, paste0("- Type S: ", type_s), paste0("- Type M: ", type_m))
}

.report_effect_precision <- function(e) {
  decision <- function(count)
    if (.report_missing(count)) "not available" else .report_rate(count, e[["n_tested"]])
  width <- e[["mean_ci_width"]]
  c(paste0("### ", .report_name(e)), "",
    paste0("- True value: ", if (.report_missing(e[["true_effect"]])) "not given"
           else .report_g(e[["true_effect"]])),
    .report_justification(e),
    paste0("- Mean interval width: ", if (.report_missing(width)) "not available"
           else sprintf("%.3g", as.numeric(width))),
    paste0("- Interval wholly outside the region (a meaningful effect): ",
           decision(e[["n_meaningful"]])),
    paste0("- Interval wholly inside the region (practical equivalence): ",
           decision(e[["n_equivalent"]])))
}

.report_section_effects <- function(rec) {
  out <- c("## Focal effects", "",
           "The rates below are those of the simulated design, each over the replicates that returned an estimate for the effect.")
  for (e in rec[["effects"]])
    out <- c(out, "", if (identical(rec[["kind"]], "precision")) .report_effect_precision(e)
                      else .report_effect_power(e, rec[["analysis"]][["joint"]]))
  out
}

.report_section_null <- function(rec) {
  z <- rec[["null"]]
  if (.report_missing(z)) return(NULL)
  precision <- identical(rec[["kind"]], "precision")
  intro <- if (precision)
    "With every focal effect set to zero, the share of the replicates whose interval fell wholly outside the region of practical equivalence:"
  else if (.report_missing(z[["alpha"]])) "With every focal effect set to zero:"
  else sprintf("With every focal effect set to zero, at alpha = %s, %s:", .report_g(z[["alpha"]]),
               .report_sidedness(z[["alternative"]], rec[["analysis"]][["joint"]]))
  count <- if (precision) "n_meaningful" else "n_significant"
  lines <- vapply(z[["effects"]], function(e)
    sprintf("- %s: %s", .report_name(e), .report_rate(e[[count]], e[["n_tested"]])),
    character(1))
  c("## Power when the null hypothesis is true", "", intro, "", lines)
}

.report_replicates <- function(n_sims, n_returned) {
  if (.report_missing(n_returned)) sprintf("%s simulated", .report_int(n_sims))
  else sprintf("%s simulated, %s returned an estimate", .report_int(n_sims),
               .report_int(n_returned))
}

.report_section_accounting <- function(rec) {
  f <- rec[["fits"]]
  out <- c("## Monte Carlo accounting", "",
           paste0("- Replicates: ", .report_replicates(f[["n_sims"]], f[["n_returned"]])))
  if (!.report_missing(f[["n_converged"]]))
    out <- c(out, sprintf("- Fits: %s converged cleanly, %s singular, %s with a warning",
                          .report_int(f[["n_converged"]]), .report_int(f[["n_singular"]]),
                          .report_int(f[["n_warning"]])))
  rungs <- rec[["analysis"]][["rungs"]]
  if (!.report_missing(rungs))
    out <- c(out, paste0("- Analysed at each rung: ", paste(vapply(rungs, function(r)
      paste(r[["name"]], .report_f3(r[["share"]])), character(1)), collapse = ", ")))
  z <- rec[["null"]]
  if (!.report_missing(z))
    out <- c(out, paste0("- Null condition: ",
                         .report_replicates(z[["n_sims"]], z[["n_returned"]])))
  note <- "Each rate's Monte Carlo standard error is sqrt(p(1 - p) / n), over the n replicates that returned an estimate, and its interval is the 95% Wilson score interval."
  if (!identical(rec[["kind"]], "precision") && !isTRUE(rec[["analysis"]][["joint"]]))
    note <- paste(note, "The standard error of Type M is the standard deviation of the exaggeration ratios over the square root of their number.")
  c(out, "", note)
}

.report_section_reproducibility <- function(rec) {
  r <- rec[["reproducibility"]]
  out <- c("## Reproducibility", "", paste0("- Software: ", r[["software"]]))
  if (.report_missing(r[["spec_file"]]))
    return(c(out, "- Specification: not given, so neither its seed nor its checksum is recorded"))
  c(out, sprintf("- Specification: %s, MD5 checksum %s", r[["spec_file"]], r[["spec_md5"]]),
    paste0("- Seed: ", if (.report_missing(r[["seed"]])) "not given"
           else .report_int(r[["seed"]])))
}

# The record as Markdown: one string, its lines joined by newlines and ending in one.
.report_render <- function(rec) {
  sections <- list(.report_section_analysis(rec), .report_section_solved(rec),
                   .report_section_effects(rec), .report_section_null(rec),
                   .report_section_accounting(rec), .report_section_reproducibility(rec))
  lines <- "# Sample-size justification"
  for (s in sections) if (length(s)) lines <- c(lines, "", s)
  paste0(paste(lines, collapse = "\n"), "\n")
}
