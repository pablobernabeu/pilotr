# Coverage of solve_curve()'s interval when the points of a curve share their random numbers, and
# when each point draws its own.
#
# solve_curve() fits the points of a curve as independent binomial counts, and the interval it
# reports is nominal only when they are. Up to and including 0.3.1, every point of a sweep reused
# the specification's seed, so replicate i was drawn from the same stream at every point and
# neighbouring points shared most of their noise. Each point now draws its own replicates. This
# script measures what the sharing did to the interval's coverage and what drawing each point
# independently restores, over enough curves for a coverage figure to mean something.
# solve_curve_calibration.R beside it measures the solver on 12 curves, which can show a solved
# size's error but not a coverage rate.
#
# A full pilotr run at this size would take hours, so the curves are emulated from the draw layout
# the package's simulator follows for the two-group Gaussian design of the calibration. With no
# random effects, row i of a replicate takes the i-th normal of that replicate's stream, and the
# first n/2 rows form the first group. The script checks that layout against simulate_design()
# before it relies on it. In the "shared" arm, each replicate's stream is the same at every point
# of a curve, so the rows of a smaller design are a prefix of a larger one's. That is what a sweep
# produced when every point reused the specification's seed. In the "independent" arm, every point
# draws afresh, as every point of a sweep now does. Each curve is then solved by the package's own
# solve_curve(), so the interval measured is the one that ships.
#
# The settings are those of solve_curve_calibration.R: three grid shapes, four effect sizes and
# three target powers at 400 replicates a point, with the grid points placed by the analytic power
# they carry. An effect-size axis is added, since a sweep over an effect leaves the design, and so
# every draw, unchanged from point to point. It has 64 subjects and the effect swept from 0.45 to
# 0.95, and the effect that gives a power of 0.8 is solved with transform = "identity".
#
# Results are written to solve_curve_coverage.txt beside this file, with LF endings. The emulation
# draws from R's own generator, seeded once, so a rerun on the same platform reproduces the file.
#
# Usage: Rscript tools/calibration/solve_curve_coverage.R [curves] [output-file] [seed]
#
# The default of 1,000 curves a setting takes several minutes. At that count the binomial standard
# error of a coverage near 0.95 is about 0.007 in each setting, and the figures the documentation
# quotes are the ones this default produces.

args <- commandArgs(trailingOnly = TRUE)
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (!nzchar(here) || is.na(here)) here <- "."
root <- normalizePath(file.path(here, "..", ".."), mustWork = FALSE)
n_curves <- if (length(args) >= 1) as.integer(args[1]) else 1000L
if (is.na(n_curves) || n_curves < 10L) stop("`curves` must be an integer of 10 or more.", call. = FALSE)
out_path <- if (length(args) >= 2) args[2] else file.path(here, "solve_curve_coverage.txt")
seed <- if (length(args) >= 3) as.integer(args[3]) else 20261013L
if (is.na(seed)) stop("`seed` must be an integer.", call. = FALSE)

# Source the whole package, as solve_curve_calibration.R and tools/parity/run_r.R do, so that the
# solve measured is the one in the source tree.
src <- file.path(root, "r", "pilotr", "R")
for (f in sort(list.files(src, pattern = "[.]R$", full.names = TRUE))) source(f)

N_SIMS <- 400L
EFFECTS <- c(0.5, 0.65, 0.8, 1.0)
TARGETS <- c(0.7, 0.8, 0.9)
GRIDS <- list(
  tight  = c(0.55, 0.68, 0.79, 0.88, 0.95),
  coarse = c(0.30, 0.55, 0.78, 0.95),
  tall   = c(0.40, 0.75, 0.92, 0.98, 0.995)
)
AXIS_N <- 64L
AXIS_GRID <- c(0.45, 0.55, 0.65, 0.75, 0.85, 0.95)
AXIS_TARGET <- 0.8
LEVEL <- 0.95
ALPHA <- 0.05

.analytic_n <- function(d, power)
  2 * stats::power.t.test(delta = d, sd = 1, sig.level = ALPHA, power = power)$n
.grid_sizes <- function(d, powers)
  pmax(6, 2 * round(vapply(powers, function(p) .analytic_n(d, p), numeric(1)) / 2))

# The layout the emulation copies, checked on the package's own simulator. At one seed, the
# residuals of a design with 72 subjects must be the first 72 residuals of the same design with
# 100, and the first half of the rows must form the first group. If either failed, the "shared" arm
# below would be emulating a sweep the package never produced.
.check_layout <- function() {
  spec <- build_spec(list(name = "layout", seed = 20260820, design_kind = "between",
    n_subject = 72, factor_name = "group", lev1 = "a", lev2 = "b", intercept = 0,
    effect = 0.5, family = "gaussian", resp_name = "score", sigma = 1))
  resid <- function(n) {
    s <- spec
    s$units$subject$n <- n
    d <- simulate_design(s)
    if (!identical(as.character(d$group), rep(c("a", "b"), each = n / 2)))
      stop(sprintf("the first half of the %d rows is not the first group.", n), call. = FALSE)
    d$score - 0.5 * ifelse(d$group == "a", -0.5, 0.5)
  }
  small <- resid(72)
  large <- resid(100)
  if (max(abs(small - large[seq_along(small)])) > 1e-12)
    stop("the rows of the smaller design are not a prefix of the larger one's.", call. = FALSE)
  invisible(TRUE)
}
.check_layout()

# Power at each point of a curve. Each element of `z` is a matrix of standard normals for one point,
# one row per replicate and one column per row of the design, of which the first n/2 form the first
# group. `shift` is the effect at each point, added to the second group. The test is the pooled
# two-sample t-test of power_design().
.rates <- function(z, sizes, shift) {
  vapply(seq_along(sizes), function(j) {
    n <- sizes[j]
    h <- n %/% 2
    a <- z[[j]][, seq_len(h), drop = FALSE]
    b <- z[[j]][, (h + 1L):n, drop = FALSE] + shift[j]
    ma <- rowMeans(a)
    mb <- rowMeans(b)
    ss <- rowSums((a - ma)^2) + rowSums((b - mb)^2)
    se <- sqrt(ss / (n - 2) * (1 / h + 1 / (n - h)))
    mean(abs((mb - ma) / se) > stats::qt(1 - ALPHA / 2, n - 2))
  }, numeric(1))
}

# Draws for one curve in either arm: one stream per replicate, reused as a prefix at every point
# in the shared arm, and drawn afresh for each point in the independent one.
.draws <- function(sizes, arm) {
  if (arm == "shared") {
    z <- matrix(stats::rnorm(N_SIMS * max(sizes)), N_SIMS)
    lapply(sizes, function(n) z[, seq_len(n), drop = FALSE])
  } else {
    lapply(sizes, function(n) matrix(stats::rnorm(N_SIMS * n), N_SIMS))
  }
}

# One solve, or NULL where solve_curve() refuses the curve. A refusal is counted, not hidden.
.solve <- function(x, rates, target, transform) {
  curve <- data.frame(x = x, power = rates, n_sims = N_SIMS)
  tryCatch(solve_curve(curve, target = target, transform = transform, level = LEVEL),
           error = function(e) NULL)
}

ARMS <- c("shared", "independent")
set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
started <- Sys.time()

# Each setting keeps, per arm, the solved value, the interval and the standard error of every
# curve the solver accepted.
results <- list()
for (gname in names(GRIDS)) {
  for (d in EFFECTS) {
    sizes <- .grid_sizes(d, GRIDS[[gname]])
    keep <- list()
    for (arm in ARMS) for (tg in TARGETS) keep[[paste(arm, tg)]] <- matrix(NA_real_, n_curves, 4)
    for (k in seq_len(n_curves)) {
      for (arm in ARMS) {
        rates <- .rates(.draws(sizes, arm), sizes, rep(d, length(sizes)))
        for (tg in TARGETS) {
          s <- .solve(sizes, rates, tg, "sqrt")
          if (!is.null(s)) keep[[paste(arm, tg)]][k, ] <- c(s$value, s$lo, s$hi, s$se)
        }
      }
    }
    for (tg in TARGETS)
      results[[length(results) + 1L]] <- list(grid = gname, d = d, target = tg,
        truth = .analytic_n(d, tg), sizes = sizes,
        shared = keep[[paste("shared", tg)]], independent = keep[[paste("independent", tg)]])
    cat(sprintf("done: %-6s effect %.2f over %s\n", gname, d, paste(sizes, collapse = " ")))
  }
}

# The effect-size axis. The design is the same at every point, so in the shared arm one stream per
# replicate serves every effect, as it did when a sweep over an effect reused the seed.
axis_truth <- stats::power.t.test(n = AXIS_N / 2, sd = 1, sig.level = ALPHA,
                                  power = AXIS_TARGET)$delta
axis <- list(shared = matrix(NA_real_, n_curves, 4), independent = matrix(NA_real_, n_curves, 4))
axis_sizes <- rep(AXIS_N, length(AXIS_GRID))
for (k in seq_len(n_curves)) {
  for (arm in ARMS) {
    z <- if (arm == "shared") rep(list(matrix(stats::rnorm(N_SIMS * AXIS_N), N_SIMS)),
                                  length(AXIS_GRID))
         else .draws(axis_sizes, "independent")
    rates <- .rates(z, axis_sizes, AXIS_GRID)
    s <- .solve(AXIS_GRID, rates, AXIS_TARGET, "identity")
    if (!is.null(s)) axis[[arm]][k, ] <- c(s$value, s$lo, s$hi, s$se)
  }
}
cat("done: effect-size axis\n")
elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))

# Summaries of one arm of one setting: coverage of the analytic value, root mean squared error,
# mean interval width, mean reported standard error and the counts behind them.
.summ <- function(m, truth) {
  ok <- !is.na(m[, 1])
  m <- m[ok, , drop = FALSE]
  list(cover = mean(m[, 2] <= truth & truth <= m[, 3]), hits = sum(m[, 2] <= truth & truth <= m[, 3]),
       solved = nrow(m), refused = sum(!ok), rmse = sqrt(mean((m[, 1] - truth)^2)),
       width = mean(m[, 3] - m[, 2]), se = mean(m[, 4]))
}

.render <- function(rows, header) {
  m <- rbind(header, do.call(rbind, rows))
  w <- apply(nchar(m), 2, max)
  cols <- lapply(seq_along(w), function(j)
    vapply(m[, j], function(s) paste0(strrep(" ", w[j] - nchar(s)), s), character(1),
           USE.NAMES = FALSE))
  sub(" +$", "", do.call(paste, c(cols, list(sep = "  "))))
}
.f <- function(x, d = 3) formatC(x, format = "f", digits = d)
# The binomial standard error of a coverage rate over independent solves, one per curve.
.mcse <- function(hits, n) sqrt((hits / n) * (1 - hits / n) / n)

# A coverage rate pooled over the effect sizes and targets of one grid shape, and its standard
# error. The three targets of a curve are solved from the same simulated rates, so their outcomes
# are correlated, and a binomial standard error over all the solves would count each curve three
# times. Each curve is instead one cluster of solves, and the pooled rate, hits over solves, takes
# the delta-method standard error of a ratio of two sums over clusters.
.pooled <- function(rs, arm) {
  h <- s <- numeric(0)
  for (rk in split(rs, vapply(rs, function(r) r$d, numeric(1)))) {
    k <- nrow(rk[[1]][[arm]])
    hit <- vapply(rk, function(r) {
      m <- r[[arm]]
      !is.na(m[, 1]) & m[, 2] <= r$truth & r$truth <= m[, 3]
    }, logical(k))
    ok <- vapply(rk, function(r) !is.na(r[[arm]][, 1]), logical(k))
    h <- c(h, rowSums(hit))
    s <- c(s, rowSums(ok))
  }
  p <- sum(h) / sum(s)
  n <- length(h)
  list(solved = sum(s), cover = p, se = sqrt(sum((h - p * s)^2) / (n * (n - 1))) / mean(s))
}

lines <- c(
  "solve_curve() coverage with shared and independent replicates across grid points",
  "",
  sprintf("Produced by tools/calibration/solve_curve_coverage.R with %d curves a setting, %d", n_curves, N_SIMS),
  sprintf("replicates a point and an emulation seed of %d, under %s.", seed, R.version.string),
  "Rerun it to regenerate this file. Nothing here is timed or dated, so a rerun that changes a",
  "line has changed a result.",
  "",
  "The curves emulate the two-group Gaussian design of solve_curve_calibration.R, whose draw",
  "layout was checked against simulate_design() before any curve was drawn. In the shared arm,",
  "each replicate's stream is reused as a prefix at every point of a curve, as when every point",
  "of a sweep reused the specification's seed. In the independent arm, every point draws afresh,",
  "as every point of a sweep now does. Each curve is solved by the package's own solve_curve(),",
  sprintf("for an interval at the %g level. At %d curves a setting, a coverage of 0.95 carries a",
          LEVEL, n_curves),
  sprintf("binomial standard error of %.3f.", sqrt(0.95 * 0.05 / n_curves)),
  "",
  sprintf("Effect sizes %s at a residual standard deviation of 1, target powers %s,",
          paste(EFFECTS, collapse = ", "), paste(TARGETS, collapse = ", ")),
  sprintf("two-tailed tests at alpha %g. Sizes are total subjects across the two groups.", ALPHA),
  "Grid shapes, as the analytic power each point carries:",
  paste0("  ", names(GRIDS), ": ", vapply(GRIDS, function(g) paste(g, collapse = " "), character(1))),
  "",
  "Coverage of the analytic size, root mean squared error and mean interval width in subjects,",
  "and the solves refused, per setting.",
  "")

rows <- lapply(results, function(r) {
  sh <- .summ(r$shared, r$truth)
  ind <- .summ(r$independent, r$truth)
  c(r$grid, .f(r$d, 2), .f(r$target, 2), .f(r$truth, 1), .f(sh$cover), .f(ind$cover),
    .f(sh$rmse, 2), .f(ind$rmse, 2), .f(sh$width, 2), .f(ind$width, 2),
    as.character(sh$refused), as.character(ind$refused))
})
lines <- c(lines, .render(rows, c("grid", "effect", "target", "analytic", "cover_sh", "cover_ind",
                                  "rmse_sh", "rmse_ind", "width_sh", "width_ind", "refused_sh",
                                  "refused_ind")), "")

lines <- c(lines, "Pooled over the effect sizes and targets of each grid shape, with the standard error",
           "of each coverage rate. The three targets of a curve are solved from the same simulated",
           "rates, so the standard error treats each curve as one cluster of three solves.", "")
pooled <- lapply(names(GRIDS), function(g) {
  rs <- Filter(function(r) r$grid == g, results)
  sh <- .pooled(rs, "shared")
  ind <- .pooled(rs, "independent")
  c(g, as.character(sh$solved), .f(sh$cover), .f(sh$se),
    as.character(ind$solved), .f(ind$cover), .f(ind$se))
})
lines <- c(lines, .render(pooled, c("grid", "solves_sh", "cover_sh", "mcse_sh", "solves_ind",
                                    "cover_ind", "mcse_ind")), "")

ind_cover <- vapply(results, function(r) .summ(r$independent, r$truth)$cover, numeric(1))
sh_cover <- vapply(results, function(r) .summ(r$shared, r$truth)$cover, numeric(1))
lines <- c(lines,
  sprintf("Across the %d settings, coverage ran from %s to %s in the shared arm and from %s",
          length(results), .f(min(sh_cover)), .f(max(sh_cover)), .f(min(ind_cover))),
  sprintf("to %s in the independent arm.", .f(max(ind_cover))), "")

ash <- .summ(axis$shared, axis_truth)
aind <- .summ(axis$independent, axis_truth)
lines <- c(lines,
  sprintf("Effect-size axis: %d subjects, effects %s, solved for a power of %g with", AXIS_N,
          paste(AXIS_GRID, collapse = " "), AXIS_TARGET),
  sprintf("transform = \"identity\". The analytic effect is %s.", .f(axis_truth, 4)), "")
lines <- c(lines, .render(list(
  c("shared", as.character(ash$solved), .f(ash$cover), .f(.mcse(ash$hits, ash$solved)),
    .f(ash$rmse, 4), .f(ash$se, 4), as.character(ash$refused)),
  c("independent", as.character(aind$solved), .f(aind$cover), .f(.mcse(aind$hits, aind$solved)),
    .f(aind$rmse, 4), .f(aind$se, 4), as.character(aind$refused))),
  c("arm", "solves", "cover", "mcse", "rmse", "mean_se", "refused")), "")
lines <- c(lines,
  "Where the reported standard error falls short of the root mean squared error, the interval",
  "claims more precision than the curve holds.")

con <- file(out_path, open = "wb")           # binary: LF endings on every platform
writeLines(lines, con, sep = "\n")
close(con)

cat("\n", paste(lines, collapse = "\n"), "\n", sep = "")
cat(sprintf("written to %s in %.0f seconds\n", out_path, elapsed))
