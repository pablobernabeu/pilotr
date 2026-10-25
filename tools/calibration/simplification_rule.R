# Power and the Type I error rate under a preregistered simplification rule, on the power
# vignette's priming design.
#
# By default, power_mixed() fits one model to every replicate, the maximal model the specification
# implies, and keeps a singular fit. An analysis plan may instead simplify the model whenever its
# fit is singular: drop the correlations, then the slopes. power_mixed() simulates such a rule when
# `formula` is given as a ladder of models. The model analysed then depends on the
# replicate's data, so the rule has its own power and its own Type I error rate. Ending on random
# intercepts alone ignores how much the effect varies over subjects and items, and a test that
# ignores a large variation rejects a true null hypothesis too often (Barr et al., 2013; Matuschek
# et al., 2017). This script measures both rates on the same replicates for three analyses: the
# maximal model, the rule as the ladder c("maximal", "zcp", "intercepts") with
# fallback = "singular", and random intercepts alone.
#
# The design is the power vignette's spec_c, 24 subjects crossed with 18 items and an effect of
# 0.06 on the log scale. Its replicates are those of
# sweep_spec(spec_c, "fixed$coefficients", design_conditions(effect = 0.06), power_mixed, ...),
# which is how the vignette runs the null condition through the rule. Each condition takes its own
# seed from the specification's, and each replicate its own seed from the condition's. A replicate
# is simulated once and analysed three ways. The rule reuses the maximal fit, and the
# intercepts-only fit when it gets that far. It fits the zero-correlation model only after a
# maximal fit that is singular or failed, so a replicate costs about 2.7 fits, against about 4.3
# for three calls of power_mixed(). The fits and the rule are pilotr's own, .fit_model() and
# .rung_passes(), sourced from this source tree. For each condition, the script then checks that
# power_mixed() returns the same counts over the first replicates, for each analysis.
#
# Beside each rejection rate, the script reports the standard deviation of the estimates over the
# mean standard error the model reported for them. A ratio above 1 means that the reported
# standard errors are too small, so the test rejects a true null hypothesis more often than alpha.
# At 200 replicates, the ratio carries a relative standard error of about 5%, where a rejection
# rate near 0.05 carries one of about 30%.
#
# Results are written to simplification_rule.txt beside this file, with LF endings. Every replicate
# takes its seed from the specification, so a rerun with the same lme4 reproduces the file.
#
# Usage: Rscript tools/calibration/simplification_rule.R [replicates] [output-file]
#
# The default of 200 replicates a condition takes about 12 minutes on one core.

args <- commandArgs(trailingOnly = TRUE)
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (!nzchar(here) || is.na(here)) here <- "."
root <- normalizePath(file.path(here, "..", ".."), mustWork = FALSE)
n_sims <- if (length(args) >= 1) as.integer(args[1]) else 200L
if (is.na(n_sims) || n_sims < 20L) stop("`replicates` must be an integer of 20 or more.", call. = FALSE)
out_path <- if (length(args) >= 2) args[2] else file.path(here, "simplification_rule.txt")
for (pkg in c("lme4", "lmerTest"))
  if (!requireNamespace(pkg, quietly = TRUE))
    stop(sprintf("this calibration needs the '%s' package.", pkg), call. = FALSE)

# Source the whole package, as the other calibration scripts do, so that the fits and the rule
# measured are the ones in the source tree.
src <- file.path(root, "r", "pilotr", "R")
for (f in sort(list.files(src, pattern = "[.]R$", full.names = TRUE))) source(f)

# The power vignette's priming design, as its chunk "priming-spec" builds it.
spec <- build_spec(list(
  name = "priming", seed = 1, design_kind = "within", include_items = TRUE,
  n_subject = 24, n_item = 18,
  factor_name = "condition", lev1 = "related", lev2 = "unrelated",
  intercept = 6, effect = 0.06,
  subj_int_sd = 0.12, subj_slope_sd = 0.04, subj_corr = 0.2,
  item_int_sd = 0.08, item_slope_sd = 0.02, item_corr = -0.1,
  family = "shifted_lognormal", resp_name = "RT", sigma = 0.3, shift = 200))
EFFECT <- spec$fixed$coefficients$effect
ALPHA <- 0.05
LADDER <- c("maximal", "zcp", "intercepts")
ANALYSES <- c("maximal", "rule", "intercepts")
N_CHECK <- 10L

# The specification at each condition, as sweep_spec() builds it: the condition merged into the
# coefficients, and the seed set to the condition's own.
conditions <- design_conditions(effect = EFFECT)
grid_seeds <- replicate_seeds(spec$seed, length(conditions))
points <- lapply(seq_along(conditions), function(i) {
  point <- spec
  point$seed <- grid_seeds[i]
  point$fixed$coefficients <- utils::modifyList(spec$fixed$coefficients, conditions[[i]])
  point
})

# The effect's estimate, standard error and p-value in a fit, NA where the fit failed, as the
# replicate loop of power_mixed() reads them.
read_fit <- function(f) {
  if (is.null(f$fit)) return(c(est = NA_real_, se = NA_real_, p = NA_real_))
  co <- summary(f$fit)$coefficients
  c(est = co["effect", "Estimate"], se = co["effect", "Std. Error"], p = co["effect", "Pr(>|t|)"])
}

# One replicate of one condition, analysed three ways.
analyse <- function(point, seed) {
  s <- point
  s$seed <- seed
  d <- model_data(point, simulate_design(s, validate = FALSE))
  fit <- function(structure) .fit_model(model_formula(point, structure), d,
                                        family = "shifted_lognormal", test = TRUE, auto = TRUE)
  fits <- list(maximal = fit("maximal"), intercepts = fit("intercepts"))
  rung <- "intercepts"
  if (.rung_passes(fits$maximal, "singular")) {
    rung <- "maximal"
  } else {
    zcp <- fit("zcp")
    if (.rung_passes(zcp, "singular")) {
      rung <- "zcp"
      fits$zcp <- zcp
    }
  }
  rule <- fits[[rung]]
  out <- lapply(list(maximal = fits$maximal, rule = rule, intercepts = fits$intercepts),
                function(f) c(read_fit(f), singular = isTRUE(f$singular)))
  list(values = out, rung = rung)
}

# What power_mixed() counts over replicates `runs`, for the check against it.
counts <- function(runs, a) {
  p <- vapply(runs, function(r) r$values[[a]][["p"]], numeric(1))
  c(n_significant = sum(p < ALPHA, na.rm = TRUE), n_returned = sum(!is.na(p)),
    n_singular = sum(vapply(runs, function(r) r$values[[a]][["singular"]] == 1, logical(1))))
}

started <- Sys.time()
results <- lapply(seq_along(points), function(i) {
  point <- points[[i]]
  seeds <- replicate_seeds(point$seed, n_sims)
  runs <- lapply(seeds, function(sd) analyse(point, sd))
  # The first replicates against power_mixed() over the same specification, whose replicate seeds
  # are a prefix of these.
  formulas <- list(maximal = NULL, rule = LADDER, intercepts = "intercepts")
  for (a in ANALYSES) {
    pm <- power_mixed(point, formula = formulas[[a]], n_sims = N_CHECK)
    want <- c(n_significant = pm$n_significant[["effect"]], n_returned = pm$n_returned,
              n_singular = pm$n_singular)
    got <- counts(runs[seq_len(N_CHECK)], a)
    if (!identical(as.numeric(want), as.numeric(got)))
      stop(sprintf("the %s analysis disagrees with power_mixed() over the first %d replicates.",
                   a, N_CHECK), call. = FALSE)
    if (a == "rule") {
      kept <- vapply(runs[seq_len(N_CHECK)], function(r) r$rung, character(1))
      if (!identical(as.numeric(pm$rung_shares),
                     as.numeric(tabulate(match(kept, LADDER), 3L) / N_CHECK)))
        stop("the rule's rungs disagree with power_mixed().", call. = FALSE)
    }
  }
  cat(sprintf("done: effect %g\n", point$fixed$coefficients$effect))
  runs
})
elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))

.f <- function(x, d = 3) formatC(x, format = "f", digits = d)
.render <- function(rows, header) {
  m <- rbind(header, do.call(rbind, rows))
  w <- apply(nchar(m), 2, max)
  cols <- lapply(seq_along(w), function(j)
    vapply(m[, j], function(s) paste0(strrep(" ", w[j] - nchar(s)), s), character(1),
           USE.NAMES = FALSE))
  sub(" +$", "", do.call(paste, c(cols, list(sep = "  "))))
}

# One row per analysis and condition. The rate is the Type I error rate where the effect is zero
# and power where it is not, with the Monte Carlo error power_mixed() reports for it.
rows <- list()
for (a in ANALYSES) {
  for (i in seq_along(points)) {
    runs <- results[[i]]
    v <- do.call(rbind, lapply(runs, function(r) r$values[[a]]))
    ok <- !is.na(v[, "p"])
    rate <- .rate_with_error(sum(v[ok, "p"] < ALPHA), sum(ok), "r")
    ratio <- stats::sd(v[ok, "est"]) / mean(v[ok, "se"])
    shares <- if (a == "rule")
      paste(.f(tabulate(match(vapply(runs, function(r) r$rung, character(1)), LADDER), 3L) /
               length(runs)), collapse = " ")
    else "-"
    rows[[length(rows) + 1L]] <- c(
      a, .f(points[[i]]$fixed$coefficients$effect, 2), as.character(sum(v[ok, "p"] < ALPHA)),
      as.character(sum(ok)), .f(rate$r), .f(rate$r_mcse), sprintf("[%s, %s]", .f(rate$r_lo),
      .f(rate$r_hi)), .f(ratio), as.character(sum(v[, "singular"] == 1)), shares)
  }
}

lines <- c(
  "Power and Type I error under a preregistered simplification rule",
  "",
  sprintf("Produced by tools/calibration/simplification_rule.R with %d replicates a condition, under", n_sims),
  sprintf("%s with lme4 %s and lmerTest %s.", R.version.string,
          as.character(utils::packageVersion("lme4")),
          as.character(utils::packageVersion("lmerTest"))),
  "Rerun it to regenerate this file. Nothing here is timed or dated, so a rerun that changes a",
  "line has changed a result.",
  "",
  "The power vignette's priming design: 24 subjects crossed with 18 items, by-subject and by-item",
  "intercepts and slopes, an effect of 0.06 on the log scale, two-sided tests at alpha 0.05 with",
  "Satterthwaite's degrees of freedom. The replicates are those of sweep_spec() over",
  "design_conditions(effect = 0.06) with power_mixed(), and each is analysed three ways:",
  "",
  "  maximal     the model the specification implies, (1 + effect | subject) + (1 + effect | item),",
  "              with singular fits kept",
  "  rule        the ladder c(\"maximal\", \"zcp\", \"intercepts\") with fallback = \"singular\"",
  "  intercepts  random intercepts alone, (1 | subject) + (1 | item)",
  "",
  sprintf("Over the first %d replicates of each condition, power_mixed() returned the same counts", N_CHECK),
  "for each analysis, and the same rungs for the rule.",
  "",
  "The rate is the Type I error rate where the effect is 0 and power where it is 0.06, with its",
  "Monte Carlo standard error and 95% Wilson interval. sd_ratio is the standard deviation of the",
  "estimates over the mean standard error reported for them, which exceeds 1 when the standard",
  "errors are too small. n_singular counts the replicates whose analysed fit was singular, and",
  "the shares are those of the rule's replicates kept at the maximal, zero-correlation and",
  "intercepts-only rungs, in that order.",
  "",
  .render(rows, c("analysis", "effect", "n_sig", "n_returned", "rate", "mcse", "ci95",
                  "sd_ratio", "n_singular", "shares")))

con <- file(out_path, open = "wb")           # binary: LF endings on every platform
writeLines(lines, con, sep = "\n")
close(con)

cat("\n", paste(lines, collapse = "\n"), "\n", sep = "")
cat(sprintf("written to %s in %.0f seconds\n", out_path, elapsed))
