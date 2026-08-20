# Bayesian design analysis, emitted as a script for the user to run.
#
# Nothing in this file fits a model, draws from a posterior, or touches Stan. The whole file is
# string templating, so it works with brms absent, which is what keeps a Bayesian design
# analysis reachable from the webR build of the no-code app, where Stan cannot run at all
# because it compiles C++ at fit time.
#
# The emitted script decides about each focal effect on two criteria at once, a Savage-Dickey
# Bayes factor and a highest-density interval against a region of practical equivalence, and it
# withholds every verdict when the sampler has not converged. The frequentist analogue of the
# interval half of that rule is precision_design(), which runs in place because lme4 needs no
# compiler.

# Fill in and check one of the settings lists (`rule`, `gate`). A missing element takes its
# default, so that rule = list(bf = 3) means the default rule with a stricter threshold. Taking
# NULL instead would leave a rule with no ROPE and no interval mass.
.da_settings <- function(x, defaults, what) {
  if (!is.list(x) || (length(x) && is.null(names(x))))
    stop("`", what, "` must be a named list", call. = FALSE)
  extra <- setdiff(names(x), names(defaults))
  if (length(extra))
    stop("unknown `", what, "` element(s) ", paste(sprintf("'%s'", extra), collapse = ", "),
         "; the recognised ones are ", paste(sprintf("'%s'", names(defaults)), collapse = ", "),
         call. = FALSE)
  out <- defaults
  for (k in names(x)) out[[k]] <- x[[k]]
  for (k in names(out)) if (!.is_scalar_number(out[[k]]))
    stop("`", what, "$", k, "` must be a single finite number", call. = FALSE)
  out
}

# Resolve `focal` to coefficient names and their true values.
#
# A character vector names the effects and leaves the true values unknown, while a named
# numeric vector supplies both. The type decides which it is, and the presence of names does
# not, so that a character vector which happens to carry names is still read as a list of effect
# names.
.da_focal <- function(focal, spec) {
  if (!length(focal) || !(is.character(focal) || is.numeric(focal)))
    stop("`focal` must be a character vector of coefficient names, or a named numeric vector ",
         "mapping coefficient names to their true values", call. = FALSE)
  if (is.numeric(focal)) {
    nms <- names(focal)
    if (is.null(nms) || anyNA(nms) || !all(nzchar(nms)))
      stop("a numeric `focal` has to be fully named, since the names are the coefficients and ",
           "the values are their true values", call. = FALSE)
    true <- as.numeric(focal)
  } else {
    nms <- as.character(focal)
    if (anyNA(nms) || !all(nzchar(nms)))
      stop("`focal` contains an empty or missing coefficient name", call. = FALSE)
    true <- rep(NA_real_, length(nms))
  }
  if (anyDuplicated(nms))
    stop("duplicated focal effect(s): ",
         paste(sprintf("'%s'", unique(nms[duplicated(nms)])), collapse = ", "), call. = FALSE)
  # A focal name that is not a fixed coefficient of this design produces a script that stops at
  # its first hypothesis test, which on a cluster is discovered hours into a queued run. Saying
  # so at emission time costs nothing.
  avail <- names(spec$fixed$coefficients)
  miss <- setdiff(nms, avail)
  if (length(miss))
    warning(sprintf(
      "focal effect%s %s %s not a fixed coefficient of this design, so the emitted script will stop at %s hypothesis test; the design's coefficients are %s",
      if (length(miss) > 1) "s" else "", paste(sprintf("'%s'", miss), collapse = ", "),
      if (length(miss) > 1) "are" else "is", if (length(miss) > 1) "those" else "that",
      if (length(avail)) paste(sprintf("'%s'", avail), collapse = ", ") else "(none)"),
      call. = FALSE)
  list(names = nms, true = stats::setNames(true, nms))
}

# A named numeric vector as an R source literal, at the full precision of .num_literal(). The
# .r_literal() walk drops the names of an atomic vector, and the true values have to stay
# attached to the coefficients they belong to.
.da_named_num <- function(x) {
  parts <- sprintf("%s = %s", vapply(names(x), .r_name, character(1)),
                   vapply(unname(x), .num_literal, character(1)))
  paste0("c(", paste(parts, collapse = ", "), ")")
}

# The formula, family and priors for the design, taken from brms_bridge().
#
# brms_bridge() prints its ready-to-fit model as a side effect, which is useful at the console
# and unwanted in the middle of a string being assembled, so the printing is diverted while the
# returned list is read. Recomputing the model here instead would give the emitted script a
# second copy of the bridge's logic, free to drift away from it.
.bridge_quietly <- function(spec) {
  con <- file(tempfile(), open = "wt")
  on.exit({ sink(); close(con) }, add = TRUE)
  sink(con)
  brms_bridge(spec)
}

# A banner separating the emitted files. The marker is a comment in both R and bash, so the
# same form works above an R part and above a shell part, and it is fixed text so that a caller
# can split the output on it.
.da_banner <- function(i, n, name, lang) {
  head <- sprintf("# ===== FILE %d of %d: %s (%s) ", i, n, name, lang)
  c("", paste0(head, strrep("=", max(3L, 79L - nchar(head)))), "")
}

# The analysis script, taking one replicate from the embedded specification through to a verdict
# per focal effect. `f` carries the focal names and true values, `bridge` the model.
.da_analysis_lines <- function(spec, f, rule, gate, bridge) {
  c(
    "#!/usr/bin/env Rscript",
    "# ---------------------------------------------------------------------------",
    sprintf("# Bayesian design analysis for the pilotr design '%s'.", spec$name),
    "#",
    "# One run is one replicate. The script simulates from the specification embedded",
    "# below, fits the confirmatory model with brms, and decides about each focal",
    "# effect on a Savage-Dickey Bayes factor together with a highest-density interval",
    "# against a region of practical equivalence. No verdict is reported unless the",
    "# sampler has converged.",
    "#",
    "# Set PILOTR_REP to a replicate index to reseed the design, and PILOTR_OUTDIR to a",
    "# directory to have the replicate's summary written there as an RDS. Leave both",
    "# unset for a single interactive run.",
    "# ---------------------------------------------------------------------------",
    "",
    "library(brms)",
    "",
    "# ---- design ----",
    "",
    "# The specification is embedded at full precision, so that the script reproduces the",
    "# design on its own, with no file to read, and a coefficient the user typed as 0.3 is",
    "# still exactly 0.3 here.",
    paste0("spec <- ", .r_literal(spec)),
    "",
    "# ---- decision rule ----",
    "",
    "# bf_threshold applies in both directions. Evidence for the effect above it, with the",
    "# interval clear of the ROPE, gives a verdict of 'supported'. Evidence for the null",
    "# above it, with the interval inside the ROPE, gives 'null'. rope is a half-width on",
    "# the scale of the model's coefficients, which for the lognormal families is the log",
    "# scale, so it bounds a proportional difference, with no milliseconds involved.",
    sprintf("bf_threshold  <- %s", .num_literal(rule$bf)),
    sprintf("rope          <- %s", .num_literal(rule$rope)),
    sprintf("ci_mass       <- %s", .num_literal(rule$ci)),
    sprintf("max_rhat      <- %s", .num_literal(gate$max_rhat)),
    sprintf("max_divergent <- %s", .num_literal(gate$max_divergent)),
    "",
    sprintf("focal      <- %s", .r_literal(f$names)),
    sprintf("focal_true <- %s   # NA where the true value was not supplied",
            .da_named_num(f$true)),
    "",
    "# ---- run identity ----",
    "",
    "# Every verdict below is a function of the rule, the gate and the specification, and a results",
    "# directory pools whatever was written into it. The identity therefore travels with every row,",
    "# so that two array runs with different ROPEs cannot be combined into one table with nothing to",
    "# tell them apart. The package version is recorded too, because the generative core changed",
    "# numerically at 0.3, so the same specification and seed do not give the same data across it.",
    "spec_fingerprint <- local({",
    "  f <- tempfile()",
    "  on.exit(unlink(f), add = TRUE)",
    '  con <- file(f, open = "wb")      # binary, so the digest does not depend on the platform',
    '  writeLines(pilotr::spec_json(spec), con, sep = "\\n")',
    "  close(con)",
    "  unname(tools::md5sum(f))",
    "})",
    'pilotr_version <- as.character(utils::packageVersion("pilotr"))',
    "",
    "# ---- replicate index ----",
    "",
    "# The replicate seeds come from pilotr's own rule. Index arithmetic on the specification's",
    "# seed will not do, because consecutive seeds are not independent streams in this generator,",
    "# so from 0.3 the package derives the whole set of replicate seeds from that seed and takes",
    "# the index-th one. Using the same rule here means replicate 7 on the cluster is the",
    "# same data set as replicate 7 in power_mixed(), and one task can be reproduced on a laptop.",
    'rep_id <- suppressWarnings(as.integer(Sys.getenv("PILOTR_REP", "1")))',
    "if (is.na(rep_id) || rep_id < 1L) rep_id <- 1L",
    "spec$seed <- pilotr::replicate_seeds(spec$seed, rep_id)[rep_id]",
    'outdir <- Sys.getenv("PILOTR_OUTDIR", "")',
    "",
    "# ---- data ----",
    "",
    "d   <- pilotr::simulate_design(spec)",
    "dat <- pilotr::model_data(spec, d)   # adds the numeric contrast and interaction columns",
    "",
    "# ---- fit ----",
    "",
    '# sample_prior = "yes" is what makes the Bayes factor computable. The Savage-Dickey ratio',
    "# is a ratio of prior to posterior density at zero, so without prior draws the",
    "# denominator cannot be read off and hypothesis() has nothing to divide by. It follows",
    "# that the prior carries the Bayes factor: widening it moves the factor towards the",
    "# null, and any reported factor is a statement about this prior as much as about the",
    "# data (Kass & Raftery, 1995).",
    "fit <- brm(",
    paste0("  ", bridge$formula, ","),
    "  data   = dat,",
    paste0("  family = ", bridge$family, ","),
    "  prior  = c(",
    paste0("    ", paste(bridge$priors, collapse = ",\n    ")),
    "  ),",
    '  sample_prior = "yes",',
    "  chains = 4, iter = 4000, warmup = 2000, cores = 4,",
    "  control = list(adapt_delta = 0.95),",
    "  seed = spec$seed",
    ")",
    "",
    "# ---- convergence gate ----",
    "",
    "# Checked before any verdict is formed. A conclusion drawn from a fit that has not",
    "# converged is worse than no conclusion, because it carries the authority of a number",
    "# without the sampling behind it. R-hat here is the rank-normalised version of Vehtari",
    "# et al. (2021), which is why the limit is near 1.01 where the older one was 1.1, and",
    "# divergent transitions are counted as a share of post-warmup draws, since a handful in",
    "# a long run means something different from a handful in a short one.",
    "draws <- brms::as_draws_df(fit)",
    'rhats <- posterior::summarise_draws(draws, "rhat")$rhat',
    "obs_rhat <- max(rhats[is.finite(rhats)])",
    "np  <- brms::nuts_params(fit)",
    'div <- np$Value[np$Parameter == "divergent__"]',
    "obs_divergent <- if (length(div)) sum(div) / length(div) else 0",
    "gate_ok <- obs_rhat <= max_rhat && obs_divergent <= max_divergent",
    "gate_msg <- paste0(",
    '  sprintf("the convergence gate failed. Maximum R-hat was %.4f against a limit of %.4f",',
    "          obs_rhat, max_rhat),",
    '  sprintf(", and %.2f%% of draws diverged against a limit of %.2f%%.",',
    "          100 * obs_divergent, 100 * max_divergent),",
    '  " No verdict is reported.")',
    "",
    "# ---- highest-density interval ----",
    "",
    "# The narrowest interval spanning ci_mass of the ordered draws. It is contiguous by",
    "# construction, so for a multimodal posterior it is the narrowest single interval rather",
    "# than the union of the dense regions. The posterior of a fixed effect in these designs",
    "# is unimodal, which is the case this estimator is for.",
    "hdi <- function(x, mass) {",
    "  x <- sort(x[is.finite(x)])",
    "  n <- length(x)",
    "  k <- max(1L, floor(mass * n))",
    "  i <- seq_len(n - k + 1L)",
    "  j <- which.min(x[i + k - 1L] - x[i])",
    "  c(lower = x[[j]], upper = x[[j + k - 1L]])",
    "}",
    "",
    "# ---- verdict per focal effect ----",
    "",
    "rows <- lapply(focal, function(effect) {",
    '  col <- paste0("b_", effect)',
    "  if (!col %in% names(draws))",
    '    stop(sprintf("\'%s\' is not a coefficient of the fitted model; its coefficients are %s",',
    '                 effect, paste(sprintf("\'%s\'", sub("^b_", "",',
    '                   grep("^b_", names(draws), value = TRUE))), collapse = ", ")),',
    "         call. = FALSE)",
    "  post <- draws[[col]]",
    "",
    "  # The Savage-Dickey density ratio (Wagenmakers et al., 2010). For a point null,",
    "  # hypothesis() reports Evid.Ratio as the posterior density at zero over the prior",
    "  # density at zero, which is evidence FOR the null. The Bayes factor for the effect is",
    "  # therefore its reciprocal. Both are kept and named, because getting the direction",
    "  # wrong inverts every conclusion below while leaving the numbers looking plausible.",
    '  h <- brms::hypothesis(fit, paste(effect, "= 0"))',
    "  bf_null   <- h$hypothesis$Evid.Ratio",
    "  bf_effect <- if (is.finite(bf_null) && bf_null > 0) 1 / bf_null else NA_real_",
    "",
    "  itv <- hdi(post, ci_mass)",
    '  outside <- itv[["lower"]] >  rope || itv[["upper"]] < -rope',
    '  inside  <- itv[["lower"]] > -rope && itv[["upper"]] <  rope',
    "",
    "  # Both criteria have to agree before the design analysis commits either way, and the",
    "  # third answer stays available when they do not (Kruschke, 2018). isTRUE() guards the",
    "  # comparisons, so a Bayes factor that came back NA leaves the verdict inconclusive",
    "  # and no missing value propagates into it.",
    "  verdict <- if (!gate_ok) NA_character_",
    '    else if (isTRUE(bf_effect > bf_threshold) && outside) "supported"',
    '    else if (isTRUE(bf_null   > bf_threshold) && inside)  "null"',
    '    else "inconclusive"',
    "",
    "  data.frame(rep = rep_id, param = effect, true = focal_true[[effect]],",
    "             estimate = stats::median(post),",
    '             lower = itv[["lower"]], upper = itv[["upper"]],',
    "             bf_effect = bf_effect, bf_null = bf_null,",
    "             max_rhat = obs_rhat, p_divergent = obs_divergent,",
    "             verdict = verdict,",
    "             bf_threshold = bf_threshold, rope = rope, ci_mass = ci_mass,",
    "             max_rhat_limit = max_rhat, max_divergent_limit = max_divergent,",
    "             spec_fingerprint = spec_fingerprint, pilotr_version = pilotr_version,",
    "             stringsAsFactors = FALSE)",
    "})",
    "",
    "# ---- summary ----",
    "",
    "summary_table <- do.call(rbind, rows)",
    "row.names(summary_table) <- NULL",
    "print(summary_table, digits = 3)",
    "if (!gate_ok) message(gate_msg)",
    "",
    "if (nzchar(outdir)) {",
    "  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)",
    '  out <- file.path(outdir, sprintf("design_analysis_rep%04d.rds", rep_id))',
    "  saveRDS(summary_table, out)",
    '  cat(sprintf("wrote %s\\n", out))',
    "}"
  )
}

# The SLURM array wrapper, one replicate per task. It has to run under any user's account
# on any cluster, so the two site-specific values, the account to charge and a writable
# project directory, are emitted as marked placeholders for the user to fill in, and the
# analysis script is invoked from wherever the parts were saved, with no fixed path anywhere
# in it. The account stays an #SBATCH directive because the scheduler reads directives
# before the shell runs, so a directive cannot take its value from a shell variable.
.da_slurm_lines <- function(spec) {
  c(
    "#!/bin/bash",
    "# =============================================================================",
    "# SLURM array job: Bayesian design analysis for the pilotr design",
    sprintf("# '%s'.", spec$name),
    "#",
    "# One array task per replicate. Each task simulates its own data set, fits the brms",
    "# model across its cores, and writes one RDS to $PROJECT_DIR/results, which",
    "# aggregate_design_analysis.R then combines.",
    "#",
    "# The R library at $PROJECT_DIR/Rlib has to carry pilotr, brms and a working Stan",
    "# backend, installed once. Nothing here compiles Stan on the fly, and a task that",
    "# has to build the model from scratch will spend most of its wall time on it.",
    "#",
    "# Submit from the directory the three parts were saved into, so that this wrapper",
    "# finds design_analysis.R beside itself:",
    "#   sbatch design_analysis.slurm                                    # 100 replicates",
    "#   sbatch --array=1 --partition=devel design_analysis.slurm        # smoke test",
    "# =============================================================================",
    "# ---- EDIT THESE two values before submitting --------------------------------",
    "# 1. The account (allocation) the scheduler charges. The scheduler reads #SBATCH",
    "#    directives before the shell runs, so this one cannot come from a variable.",
    "#SBATCH --account=EDIT_ME_ACCOUNT",
    "# 2. PROJECT_DIR, just after the directives below: a directory you can write to,",
    "#    holding the R library at Rlib/.",
    "# -----------------------------------------------------------------------------",
    "#SBATCH --job-name=pilotr_bda",
    "#SBATCH --array=1-100               # one replicate per task",
    "#SBATCH --partition=short           # short <=12h",
    "#SBATCH --nodes=1",
    "#SBATCH --ntasks=1",
    "#SBATCH --cpus-per-task=4           # one core per brms chain",
    "#SBATCH --mem=16G",
    "#SBATCH --time=04:00:00",
    "",
    "set -uo pipefail",
    "conda deactivate 2>/dev/null || true",
    "module purge 2>/dev/null || true",
    "module load R    # or the versioned R module your cluster names",
    "",
    "PROJECT_DIR=$HOME/pilotr_toolkit    # EDIT THIS: writable, holds Rlib/ (see above)",
    "",
    "# The analysis script is the one saved next to this wrapper. sbatch runs a spooled",
    "# copy of the wrapper, so $0 does not name the saved location; SLURM_SUBMIT_DIR is",
    "# the directory sbatch was invoked from, which the submit instructions above pin to",
    "# the directory holding the parts. $0 covers running the wrapper directly.",
    'SCRIPT_DIR=${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}',
    "",
    'export R_LIBS=$PROJECT_DIR/Rlib:${R_LIBS:-}   # pilotr, brms and the Stan backend',
    "export PILOTR_OUTDIR=$PROJECT_DIR/results/design_analysis",
    "export PILOTR_REP=$SLURM_ARRAY_TASK_ID",
    "export OMP_NUM_THREADS=1                      # one thread per chain, the chains parallelise",
    'mkdir -p "$PILOTR_OUTDIR" "$PROJECT_DIR/tmp"',
    "export TMPDIR=$PROJECT_DIR/tmp",
    "",
    'echo "Host $(hostname) | replicate $SLURM_ARRAY_TASK_ID | cpus ${SLURM_CPUS_PER_TASK} | $(date)"',
    'Rscript "$SCRIPT_DIR/design_analysis.R"',
    'echo "Exit $? | $(date)"'
  )
}

# The aggregator, which reads the per-task RDS files and rebuilds the same summary table across
# replicates.
.da_aggregate_lines <- function() {
  c(
    "#!/usr/bin/env Rscript",
    "# ---------------------------------------------------------------------------",
    "# Combine the per-task RDS files written by design_analysis.R into one table.",
    "# Env: PILOTR_OUTDIR, the directory the array tasks wrote to.",
    "# ---------------------------------------------------------------------------",
    "",
    'outdir <- Sys.getenv("PILOTR_OUTDIR", ".")',
    'files <- list.files(outdir, pattern = "^design_analysis_rep[0-9]+[.]rds$",',
    "                    full.names = TRUE)",
    'if (!length(files)) stop("no per-replicate RDS files under ", outdir, call. = FALSE)',
    "",
    "# A results directory is a pool of whatever was written into it, and every verdict is a",
    "# function of the rule, the gate and the specification. Combining two array runs that differed",
    "# in any of those produces one table with nothing to tell them apart, so the run identity each",
    "# replicate carries is checked before anything is bound together. Give each run its own",
    "# directory; pooling across rules is a decision to make deliberately, not by globbing.",
    "parts <- lapply(files, readRDS)",
    "if (length(unique(lapply(parts, names))) > 1L)",
    '  stop("the per-replicate files under ", outdir, " do not share a column set, so they were ",',
    '       "written by different versions of this script. Aggregate each run on its own.",',
    "       call. = FALSE)",
    "reps <- do.call(rbind, parts)",
    "",
    'run_cols <- c("bf_threshold", "rope", "ci_mass", "max_rhat_limit", "max_divergent_limit",',
    '              "spec_fingerprint", "pilotr_version")',
    "mixed <- run_cols[vapply(run_cols, function(k) length(unique(reps[[k]])) > 1L, logical(1))]",
    "if (length(mixed))",
    '  stop("the per-replicate files under ", outdir, " come from more than one run: ",',
    '       paste(vapply(mixed, function(k) sprintf("%s takes the values %s", k,',
    '         paste(sort(unique(as.character(reps[[k]]))), collapse = ", ")), character(1)),',
    '         collapse = "; "),',
    '       ". Aggregate each run on its own.", call. = FALSE)',
    "",
    'write.csv(reps, file.path(outdir, "design_analysis_replicates.csv"), row.names = FALSE)',
    "",
    "# One row per focal effect, carrying the columns of a single replicate's table plus the",
    "# share of replicates that reached each verdict. The averaged bounds describe the",
    "# interval this design buys, and the shares are the design analysis itself, since a rule",
    "# that decides in a small share of replicates is a rule the design cannot afford.",
    "# Replicates that failed the convergence gate are counted on their own rather than",
    "# folded into 'inconclusive', because a gate failure says something about the sampler",
    "# and not about the evidence.",
    "share <- function(v, what) mean(!is.na(v) & v == what)",
    "by_param <- split(reps, factor(reps$param, levels = unique(reps$param)))",
    "agg <- do.call(rbind, lapply(by_param, function(r) data.frame(",
    "  param = r$param[1], true = r$true[1], n_reps = nrow(r),",
    "  n_gate_failed = sum(is.na(r$verdict)),",
    "  estimate = mean(r$estimate), lower = mean(r$lower), upper = mean(r$upper),",
    "  bf_effect = stats::median(r$bf_effect, na.rm = TRUE),",
    '  p_supported    = share(r$verdict, "supported"),',
    '  p_null         = share(r$verdict, "null"),',
    '  p_inconclusive = share(r$verdict, "inconclusive"),',
    "  # Constant across the rows by the check above, so the summary says which rule produced it.",
    "  bf_threshold = r$bf_threshold[1], rope = r$rope[1], ci_mass = r$ci_mass[1],",
    "  max_rhat_limit = r$max_rhat_limit[1], max_divergent_limit = r$max_divergent_limit[1],",
    "  spec_fingerprint = r$spec_fingerprint[1], pilotr_version = r$pilotr_version[1],",
    "  stringsAsFactors = FALSE)))",
    "row.names(agg) <- NULL",
    "print(agg, digits = 3)",
    'write.csv(agg, file.path(outdir, "design_analysis_summary.csv"), row.names = FALSE)',
    'cat(sprintf("combined %d replicate file(s) from %s\\n", length(files), outdir))'
  )
}

#' Generate a Bayesian design-analysis script from a specification
#'
#' Emit a runnable R script that simulates from a design specification, fits the confirmatory
#' Bayesian model that [brms_bridge()] derives for it, and decides about each focal effect on a
#' Savage-Dickey Bayes factor together with a highest-density interval tested against a region
#' of practical equivalence. The script is returned as a character string and nothing is fitted
#' here, so the function needs neither `brms` nor Stan installed.
#'
#' @details
#' The function emits a script for the user to run, for two reasons that both come down to
#' where a Stan model can be built. The no-code application ships as a webR build that
#' runs entirely in the browser, and Stan compiles C++ at fit time, so a Bayesian fit cannot run
#' there at all. A function that produced its analysis only when `brms` was present would then be
#' missing from the build most users meet first. The second reason is maintenance. Depending on
#' `brms` would pull a Stan toolchain into this package's own test and check matrix, and one
#' maintainer cannot keep that working across the platforms CRAN builds on. Emitting a script
#' keeps the analysis reproducible and open to inspection while leaving the fit where a compiler
#' is available.
#'
#' The verdict rests on two criteria that answer different questions. A Bayes factor compares the
#' null with the alternative and so reports which of the two the data favour (Kass and Raftery,
#' 1995), computed here as the Savage-Dickey density ratio, the ratio of prior to posterior
#' density at zero (Wagenmakers et al., 2010). The interval against a region of practical
#' equivalence asks instead whether the effect is large enough to matter (Kruschke, 2018).
#' Requiring the two to agree makes `"supported"` and `"null"` harder to reach than either
#' criterion alone would, and it leaves the third answer of `"inconclusive"` available when they
#' disagree, which is the answer a small pilot most often deserves.
#'
#' Because the Savage-Dickey ratio is read off the prior as well as the posterior, the emitted
#' `brm()` call sets `sample_prior = "yes"`, and the Bayes factor it produces is a statement
#' about the prior that [brms_bridge()] supplies as much as about the data. Widening that prior
#' moves the factor towards the null.
#'
#' The convergence gate is checked before any verdict is formed. R-hat is the rank-normalised
#' version of Vehtari et al. (2021), for which a limit near 1.01 is appropriate where the older
#' one was 1.1, and divergent transitions are counted as a share of post-warmup draws so that the
#' threshold means the same thing whatever the run length. When either limit is exceeded the
#' emitted script reports `NA` for every verdict and prints why, since a conclusion from a fit
#' that has not converged carries the authority of a number without the sampling behind it.
#'
#' A verdict is a function of the whole rule, the whole gate and the specification, so each
#' replicate's record carries all of them: `bf_threshold`, `rope`, `ci_mass`, `max_rhat_limit`,
#' `max_divergent_limit`, an MD5 fingerprint of the canonical specification JSON, and the pilotr
#' version that produced the data. A results directory collects whatever was written into it, and
#' the emitted aggregator globs it, so without those columns two array runs under different
#' regions of practical equivalence would combine into one table with nothing to tell them apart.
#' The aggregator stops when they disagree, and pools nothing.
#'
#' @param spec A design specification (path or list).
#' @param focal A character vector of focal coefficient names, or a named numeric vector mapping
#'   those names to their true values. The names are the coefficients of the emitted model, so an
#'   interaction is written `a:b` as in the specification, not `a_b` as in the `lme4` formula
#'   from [model_formula()]. A name that is not a fixed coefficient of the design is reported as
#'   a warning here, well before the emitted script would fail on it.
#' @param rule The decision rule, a named list with elements `bf` (the Bayes-factor threshold,
#'   applied in both directions), `rope` (the half-width of the region of practical equivalence,
#'   on the scale of the model's coefficients) and `ci` (the mass of the highest-density
#'   interval). A missing element takes its default. Setting `rope` to 0 leaves a verdict of
#'   `"null"` unreachable, since no interval of positive width lies inside a region of no width,
#'   which turns the rule into a Bayes factor plus a sign requirement.
#' @param engine The fitting engine for the emitted script. Only `"brms"` is supported.
#' @param gate The convergence gate, a named list with elements `max_rhat` (the largest
#'   acceptable R-hat over all parameters) and `max_divergent` (the largest acceptable share of
#'   post-warmup draws ending in a divergent transition). A missing element takes its default.
#' @param array Either `"none"` for the analysis script alone, or `"slurm"` to append a SLURM
#'   array wrapper that runs one replicate per task and an aggregator that combines the
#'   per-task results. The three parts are separated by `# ===== FILE n of 3` banners and are
#'   meant to be split into three files, since the middle part is shell and so is not valid R.
#'   The wrapper is written for a generic SLURM cluster and carries two placeholders, marked
#'   `EDIT` in its header, that must be filled in before submission: the `#SBATCH --account`
#'   directive and the writable `PROJECT_DIR`. It runs the `design_analysis.R` saved next to
#'   it, so submit from the directory the parts were saved into.
#' @param file Optional path. When given, the script is also written there, byte for byte
#'   with LF line endings and a single trailing newline on every platform, and returned
#'   invisibly.
#' @return A length-one character string holding the emitted script, invisibly when `file` is
#'   given. For `array = "slurm"` the string holds three banner-separated parts, an R analysis
#'   script, a bash array wrapper, and an R aggregator; the wrapper is not submittable as
#'   emitted, since its `--account` and `PROJECT_DIR` placeholders must be filled in first.
#' @references Kass, R. E. and Raftery, A. E. (1995). Bayes factors. \emph{Journal of the
#'   American Statistical Association}, 90(430), 773-795. \doi{10.1080/01621459.1995.10476572}
#'
#'   Kruschke, J. K. (2018). Rejecting or accepting parameter values in Bayesian estimation.
#'   \emph{Advances in Methods and Practices in Psychological Science}, 1(2), 270-280.
#'   \doi{10.1177/2515245918771304}
#'
#'   Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and Burkner, P.-C. (2021).
#'   Rank-normalization, folding, and localization: An improved R-hat for assessing convergence
#'   of MCMC. \emph{Bayesian Analysis}, 16(2), 667-718. \doi{10.1214/20-BA1221}
#'
#'   Wagenmakers, E.-J., Lodewyckx, T., Kuriyal, H. and Grasman, R. (2010). Bayesian hypothesis
#'   testing for psychologists: A tutorial on the Savage-Dickey method. \emph{Cognitive
#'   Psychology}, 60(3), 158-189. \doi{10.1016/j.cogpsych.2009.12.001}
#' @seealso [brms_bridge()] for the model the script fits, [precision_design()] for the
#'   frequentist analogue of the interval criterion, which runs in place, and
#'   [generate_r_script()] for the simulation-only script.
#' @examples
#' # Emitting the script needs neither brms nor Stan, which is what lets it work in the
#' # browser build of the no-code app.
#' spec   <- load_spec(pilotr_example("crossed_mixed_rt"))
#' script <- generate_design_analysis(spec, focal = c(cond = 0.05))
#' cat(head(strsplit(script, "\n")[[1]], 15), sep = "\n")
#'
#' # A stricter rule, with a wrapper for a SLURM array and its aggregator.
#' cluster <- generate_design_analysis(spec, focal = c(cond = 0.05),
#'                                     rule = list(bf = 30, rope = 0.02),
#'                                     array = "slurm")
#' cat(grep("^# ===== FILE", strsplit(cluster, "\n")[[1]], value = TRUE), sep = "\n")
#' @export
generate_design_analysis <- function(spec, focal,
                                     rule   = list(bf = 10, rope = 0.05, ci = 0.95),
                                     engine = "brms",
                                     gate   = list(max_rhat = 1.01, max_divergent = 0.01),
                                     array  = c("none", "slurm"),
                                     file   = NULL) {
  spec   <- .as_spec(spec)
  engine <- match.arg(engine, "brms")
  array  <- match.arg(array)
  rule   <- .da_settings(rule, list(bf = 10, rope = 0.05, ci = 0.95), "rule")
  gate   <- .da_settings(gate, list(max_rhat = 1.01, max_divergent = 0.01), "gate")

  if (rule$bf <= 1)
    stop("`rule$bf` must be greater than 1, since a threshold of 1 or less is met by evidence ",
         "pointing either way", call. = FALSE)
  if (rule$rope < 0) stop("`rule$rope` is a half-width, so it cannot be negative", call. = FALSE)
  if (rule$ci <= 0 || rule$ci >= 1)
    stop("`rule$ci` is the mass of an interval, so it must lie between 0 and 1", call. = FALSE)
  if (gate$max_rhat < 1)
    stop("`gate$max_rhat` cannot be below 1, which is the value R-hat approaches from above ",
         "as the chains mix", call. = FALSE)
  if (gate$max_divergent < 0 || gate$max_divergent > 1)
    stop("`gate$max_divergent` is a share of post-warmup draws, so it must lie between 0 and 1",
         call. = FALSE)

  f      <- .da_focal(focal, spec)
  bridge <- .bridge_quietly(spec)
  analysis <- .da_analysis_lines(spec, f, rule, gate, bridge)

  lines <- if (identical(array, "none")) analysis else c(
    "# =============================================================================",
    sprintf("# Bayesian design analysis for the pilotr design '%s', as three files for a",
            spec$name),
    "# SLURM array. Split this output at the 'FILE n of 3' banners below and save each",
    "# part under the name its banner gives. The first and third parts are R and run as",
    "# they stand. The middle part is bash, so it is not valid R and has to be separated",
    "# out before either part is used.",
    "# =============================================================================",
    .da_banner(1L, 3L, "design_analysis.R", "R"), analysis,
    .da_banner(2L, 3L, "design_analysis.slurm", "bash"), .da_slurm_lines(spec),
    .da_banner(3L, 3L, "aggregate_design_analysis.R", "R"), .da_aggregate_lines())

  out <- paste0(paste(lines, collapse = "\n"), "\n")
  if (!is.null(file)) {
    # Binary, so the script reaches the disk exactly as returned: LF line endings and a
    # single trailing newline on every platform. Text-mode writeLines() writes CRLF on
    # Windows, which turns the slurm part's first line into "#!/bin/bash\r", a shebang no
    # cluster can execute, and it appends a newline to a string that already ends in one.
    con <- file(file, open = "wb")
    on.exit(close(con), add = TRUE)
    writeLines(out, con, sep = "")
    return(invisible(out))
  }
  out
}
