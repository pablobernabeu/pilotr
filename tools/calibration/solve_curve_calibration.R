# Calibration of solve_curve() against an analytic answer.
#
# solve_curve() defends two design choices with measured figures: a probit link rather than a
# logit, and the delta-method interval rather than Fieller's. This script is where those figures
# come from, so that nothing quoted in the documentation rests on a run nobody can repeat. It
# builds power curves for two-group Gaussian designs, the one case stats::power.t.test() answers
# in closed form, solves each curve with the package's own solve_curve(), and compares the solved
# sample size with the analytic one.
#
# The experiment is three grid shapes by four effect sizes by three target powers: 36 checks over
# 12 curves. Grid points are placed by the analytic power they carry rather than at fixed sample
# sizes, because the size a design needs varies by a factor of four across the effect sizes swept,
# and one fixed grid would straddle the target for some of them and extrapolate for the rest. The
# three shapes differ in what they span. "tight" sits close around the targets, "coarse" is four
# widely spaced points, and "tall" reaches a power of 0.995, which is the region where a probit and
# a logit disagree most, since only the probit is linear there in the square root of the size.
#
# The curves come from sweep_spec() over power_design(), the path a user takes, so the fit sees
# exactly the data frame the package produces, including the replicate counts that weight it. The
# probit solve is solve_curve() itself and the probit fit behind Fieller's interval is the
# package's own .solve_irls(), so the comparison measures what ships rather than a copy of it. The
# logit arm is written out here, since the package has no logit to call: it is the same iteratively
# reweighted fit, the same delta method and the same heterogeneity factor, differing only in the
# link, which is what makes the two arms comparable.
#
# Every specification carries its own seed, so a rerun on the same machine reproduces the artefact
# byte for byte. Results are written to solve_curve_calibration.txt beside this file, with LF
# endings, so the figures can be read without rerunning anything.
#
# Usage: Rscript tools/calibration/solve_curve_calibration.R [replicates] [output-file] [seed]
#
# The default of 400 replicates a point runs in about a minute, nearly all of it in the
# simulation. The figures quoted in solve_curve.R and in the changelogs are the ones this default
# produces, and they are meaningless without it: the error of a solved size is dominated by Monte
# Carlo noise and falls with the square root of the replicate count. A larger count is a different
# experiment, so raise it to study the solver, not to restate the documentation.

args <- commandArgs(trailingOnly = TRUE)
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (!nzchar(here) || is.na(here)) here <- "."
root <- normalizePath(file.path(here, "..", ".."), mustWork = FALSE)
n_sims <- if (length(args) >= 1) as.integer(args[1]) else 400L
if (is.na(n_sims) || n_sims < 20L) stop("`replicates` must be an integer of 20 or more.", call. = FALSE)
out_path <- if (length(args) >= 2) args[2] else file.path(here, "solve_curve_calibration.txt")
# The seed base is an argument so that a reader can ask how much of a margin between two arms is
# the arms and how much is one draw of 12 curves.
seed_base <- if (length(args) >= 3) as.integer(args[3]) else 20260819L
if (is.na(seed_base)) stop("`seed` must be an integer.", call. = FALSE)

# Source the whole package rather than a hand-listed subset, so a new file cannot leave this
# harness measuring a stale definition. The same rule as tools/parity/run_r.R.
src <- file.path(root, "r", "pilotr", "R")
for (f in sort(list.files(src, pattern = "[.]R$", full.names = TRUE))) source(f)

EFFECTS <- c(0.5, 0.65, 0.8, 1.0)
TARGETS <- c(0.7, 0.8, 0.9)
# Grid shapes, given as the analytic power each point should carry. Every shape spans all three
# targets with a margin at both ends, because a curve whose simulated rates do not straddle the
# target is refused rather than extrapolated, and a run that measured refusals would measure the
# grid rather than the link.
GRIDS <- list(
  tight  = c(0.55, 0.68, 0.79, 0.88, 0.95),
  coarse = c(0.30, 0.55, 0.78, 0.95),
  tall   = c(0.40, 0.75, 0.92, 0.98, 0.995)
)
LEVEL <- 0.95
ALPHA <- 0.05

# Total subjects across the two groups, which is what the sweep varies and what solve_curve()
# returns. power.t.test() reports the size of one group.
.analytic_n <- function(d, power)
  2 * stats::power.t.test(delta = d, sd = 1, sig.level = ALPHA, power = power)$n

# Grid sizes are rounded to an even number so the two groups stay balanced, which is the design
# power.t.test() describes.
.grid_sizes <- function(d, powers)
  pmax(6, 2 * round(vapply(powers, function(p) .analytic_n(d, p), numeric(1)) / 2))

# The logit arm. Iteratively reweighted least squares for the binomial logit model with prior
# weights `m`, mirroring .solve_irls() line for line apart from the link, so that any difference
# between the arms is the link and not the fitting.
.logit_irls <- function(u, y, m) {
  p0 <- (m * y + 0.5) / (m + 1)
  eta <- log(p0 / (1 - p0))
  fit <- NULL; prev <- NULL
  for (it in seq_len(.SOLVE_MAXIT)) {
    mu <- 1 / (1 + exp(-eta))
    v <- mu * (1 - mu)
    w <- ifelse(v > 0, m * v, 0)
    z <- ifelse(v > 0, eta + (y - mu) / v, eta)
    fit <- .solve_wls(u, z, w)
    if (is.null(fit)) return(NULL)
    b <- c(fit$b0, fit$b1)
    eta <- fit$b0 + fit$b1 * u
    if (!is.null(prev) && max(abs(b - prev)) <= .SOLVE_TOL * (1 + max(abs(b)))) {
      fit$eta <- eta
      return(fit)
    }
    prev <- b
  }
  NULL
}

# Pearson's chi-square over the residual degrees of freedom, reported unfloored. The package floors
# the same quantity at one before it scales an interval, since a factor below one would claim more
# precision than the replicates support, but a floored value cannot show how badly a link fits.
.chisq_df <- function(mu, y, m) {
  df <- length(y) - 2L
  if (df < 1L) return(NA_real_)
  chi <- 0
  for (i in seq_along(y)) {
    v <- mu[i] * (1 - mu[i])
    if (v > 0) chi <- chi + m[i] * (y[i] - mu[i]) * (y[i] - mu[i]) / v
  }
  chi / df
}

# Delta-method solve on a fitted pair of coefficients, the inversion solve_curve() performs. Used
# for the logit arm only; the probit arm calls solve_curve() itself.
.delta_solve <- function(fit, et, z, disp) {
  us <- (et - fit$b0) / fit$b1
  g0 <- -1 / fit$b1
  g1 <- -us / fit$b1
  se <- sqrt(disp * (g0 * g0 * fit$v00 + 2 * g0 * g1 * fit$v01 + g1 * g1 * fit$v11))
  list(value = .solve_back(us, "sqrt"), lo = .solve_back(us - z * se, "sqrt"),
       hi = .solve_back(us + z * se, "sqrt"))
}

# Fieller's interval for the ratio (et - b0) / b1: the set of points whose fitted value differs
# from `et` by no more than its own standard error times z. It is bounded only where the leading
# coefficient stays positive, which is the same condition as the slope refusal in solve_curve().
.fieller <- function(fit, et, z, disp) {
  v00 <- disp * fit$v00; v01 <- disp * fit$v01; v11 <- disp * fit$v11
  a <- fit$b1 * fit$b1 - z * z * v11
  b <- fit$b1 * (fit$b0 - et) - z * z * v01
  cc <- (fit$b0 - et) * (fit$b0 - et) - z * z * v00
  disc <- b * b - a * cc
  if (a <= 0 || disc <= 0)
    return(list(lo = NA_real_, hi = NA_real_, u_lo = NA_real_, u_hi = NA_real_))
  r <- sqrt(disc)
  u_lo <- (-b - r) / a; u_hi <- (-b + r) / a
  list(lo = .solve_back(u_lo, "sqrt"), hi = .solve_back(u_hi, "sqrt"),
       u_lo = u_lo, u_hi = u_hi)
}

# The logit arm and Fieller's interval are written here rather than taken from the package, so
# they are checked before anything is measured with them. A comparison arm that is quietly wrong
# would make the shipped choice look good for the wrong reason.
.self_check <- function() {
  u <- sqrt(c(20, 40, 60, 80, 100))
  m <- rep(200, 5)
  s <- c(40, 90, 124, 156, 176)
  y <- s / m
  for (link in c("probit", "logit")) {
    fit <- if (link == "probit") .solve_irls(u, y, m) else .logit_irls(u, y, m)
    ref <- stats::glm(cbind(s, m - s) ~ u, family = stats::binomial(link = link),
                      control = stats::glm.control(epsilon = 1e-12, maxit = 100))
    b <- unname(stats::coef(ref))
    cv <- summary(ref)$cov.unscaled
    got <- c(fit$b0, fit$b1, fit$v00, fit$v01, fit$v11)
    want <- c(b[1], b[2], cv[1, 1], cv[1, 2], cv[2, 2])
    if (max(abs(got - want) / (1 + abs(want))) > 1e-7)
      stop(sprintf("the %s fit does not reproduce glm(): %s against %s", link,
                   paste(signif(got, 8), collapse = " "), paste(signif(want, 8), collapse = " ")),
           call. = FALSE)
  }
  # Fieller's bounds are the two points at which the fitted value differs from the target by
  # exactly z standard errors, so both must solve that equality.
  fit <- .solve_irls(u, y, m)
  et <- as241(0.8); zz <- as241(0.975)
  fi <- .fieller(fit, et, zz, 1)
  for (r in c(fi$u_lo, fi$u_hi)) {
    lhs <- (fit$b0 + fit$b1 * r - et)^2
    rhs <- zz * zz * (fit$v00 + 2 * r * fit$v01 + r * r * fit$v11)
    if (abs(lhs - rhs) > 1e-10 * (1 + abs(rhs)))
      stop("Fieller's bounds do not solve their own defining equation.", call. = FALSE)
  }
}
.self_check()

# Column-aligned rendering, since the artefact exists to be read rather than parsed.
.render <- function(rows, header) {
  m <- rbind(header, do.call(rbind, rows))
  w <- apply(nchar(m), 2, max)
  cols <- lapply(seq_along(w), function(j)
    vapply(m[, j], function(s) paste0(strrep(" ", w[j] - nchar(s)), s), character(1),
           USE.NAMES = FALSE))
  sub(" +$", "", do.call(paste, c(cols, list(sep = "  "))))
}

.pct <- function(x) sprintf("%.2f", x)
.num <- function(x, d = 1) formatC(x, format = "f", digits = d)

started <- Sys.time()
curves <- list()
checks <- list()
refusals <- character(0)
z <- as241((1 + LEVEL) / 2)
seed <- seed_base

for (gname in names(GRIDS)) {
  for (d in EFFECTS) {
    sizes <- .grid_sizes(d, GRIDS[[gname]])
    seed <- seed + 1L
    spec <- build_spec(list(name = "calib", seed = seed, design_kind = "between",
      n_subject = sizes[1], factor_name = "group", lev1 = "a", lev2 = "b",
      intercept = 0, effect = d, family = "gaussian", resp_name = "score", sigma = 1))
    # sweep_spec() carries the specification's seed to every grid point, so replicate i is drawn
    # from the same seed at every size and the points of one curve share their random numbers.
    # That is what a user's sweep does, and it is left alone: a calibration that broke the
    # correlation would be measuring a curve the package does not produce.
    curve <- sweep_spec(spec, "units$subject$n", sizes, power_design,
                        n_sims = n_sims, alpha = ALPHA)

    u <- sqrt(curve$n); yv <- curve$power; mv <- curve$n_sims
    pfit <- .solve_irls(u, yv, mv)
    lfit <- .logit_irls(u, yv, mv)
    if (is.null(pfit) || is.null(lfit))
      stop(sprintf("a fit did not settle on grid '%s' at an effect of %g.", gname, d), call. = FALSE)
    pchi <- .chisq_df(stats::pnorm(pfit$b0 + pfit$b1 * u), yv, mv)
    lchi <- .chisq_df(1 / (1 + exp(-lfit$eta)), yv, mv)
    curves[[length(curves) + 1L]] <- list(grid = gname, d = d, points = length(sizes),
      lo = min(sizes), hi = max(sizes), rates = yv, pchi = pchi, lchi = lchi)

    for (tg in TARGETS) {
      truth <- .analytic_n(d, tg)
      pr <- tryCatch(solve_curve(curve, target = tg, level = LEVEL),
                     error = function(e) conditionMessage(e))
      # A refusal is a result, not a crash. It is recorded and left out of the summaries, whose
      # check counts then say how many solves they rest on. The shipped replicate count produces
      # none; a much smaller one produces several, which is itself the reason not to lower it.
      if (is.character(pr)) {
        refusals <- c(refusals, sprintf("  %s, effect %.2f, target %.2f: %s", gname, d, tg, pr))
        next
      }
      lg <- .delta_solve(lfit, log(tg / (1 - tg)), z, max(lchi, 1))
      fi <- .fieller(pfit, as241(tg), z, pr$dispersion)
      checks[[length(checks) + 1L]] <- list(
        grid = gname, d = d, target = tg, truth = truth,
        p_value = pr$value, p_lo = pr$lo, p_hi = pr$hi,
        l_value = lg$value, l_lo = lg$lo, l_hi = lg$hi,
        f_lo = fi$lo, f_hi = fi$hi,
        p_err = 100 * abs(pr$value - truth) / truth,
        l_err = 100 * abs(lg$value - truth) / truth,
        p_cov = truth >= pr$lo && truth <= pr$hi,
        l_cov = truth >= lg$lo && truth <= lg$hi,
        f_cov = is.finite(fi$lo) && truth >= fi$lo && truth <= fi$hi,
        p_width = pr$hi - pr$lo,
        f_width = if (is.finite(fi$lo)) fi$hi - fi$lo else NA_real_,
        # solve_curve() refuses a point outside the swept range, so a logit that lands there is a
        # solve the package would not have reported at all.
        l_out = lg$value < min(sizes) || lg$value > max(sizes))
    }
    cat(sprintf("solved: %-6s effect %.2f over %s\n", gname, d, paste(sizes, collapse = " ")))
  }
}

.col <- function(f, mode) vapply(checks, function(r) r[[f]], mode)
num <- function(f) .col(f, numeric(1))
flag <- function(f) .col(f, logical(1))
grid_of <- .col("grid", character(1))
elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))

lines <- c(
  "solve_curve() calibration against stats::power.t.test()",
  "",
  sprintf("Produced by tools/calibration/solve_curve_calibration.R with %d replicates a point and a", n_sims),
  sprintf("seed base of %d, under %s. Rerun it to regenerate this", seed_base, R.version.string),
  "file. Nothing here is timed or dated, so a rerun that changes a line has changed a result.",
  "",
  sprintf("Effect sizes %s at a residual standard deviation of 1, target powers %s,",
          paste(EFFECTS, collapse = ", "), paste(TARGETS, collapse = ", ")),
  sprintf("two-tailed tests at alpha %g, intervals at the %g level. Sizes are total subjects", ALPHA, LEVEL),
  "across the two groups, which is twice what power.t.test() reports.",
  "",
  "Grid shapes, as the analytic power each point carries:",
  paste0("  ", names(GRIDS), ": ", vapply(GRIDS, function(g) paste(g, collapse = " "), character(1))),
  "",
  "Solved size against the analytic size, per check.",
  "")

rows <- lapply(checks, function(r) c(r$grid, .num(r$d, 2), .num(r$target, 2), .num(r$truth),
  .num(r$p_value), .pct(r$p_err), .num(r$l_value), .pct(r$l_err)))
lines <- c(lines, .render(rows, c("grid", "effect", "target", "analytic", "probit",
                                  "err%", "logit", "err%")), "")

lines <- c(lines, "Interval bounds and coverage of the analytic size, per check.", "")
rows <- lapply(checks, function(r) c(r$grid, .num(r$d, 2), .num(r$target, 2), .num(r$truth),
  .num(r$p_lo), .num(r$p_hi), if (r$p_cov) "yes" else "no",
  if (is.finite(r$f_lo)) .num(r$f_lo) else "none", if (is.finite(r$f_hi)) .num(r$f_hi) else "none",
  if (r$f_cov) "yes" else "no",
  .num(r$l_lo), .num(r$l_hi), if (r$l_cov) "yes" else "no"))
lines <- c(lines, .render(rows, c("grid", "effect", "target", "analytic", "delta_lo", "delta_hi",
                                  "cov", "fieller_lo", "fieller_hi", "cov", "logit_lo",
                                  "logit_hi", "cov")), "")

lines <- c(lines, "Fit of each link to each curve: Pearson's chi-square over its degrees of",
           "freedom, unfloored, and the simulated rates the fit saw.", "")
rows <- lapply(curves, function(r) c(r$grid, .num(r$d, 2), as.character(r$points),
  paste0(r$lo, "-", r$hi), .num(r$pchi, 2), .num(r$lchi, 2),
  paste(formatC(r$rates, format = "f", digits = 3), collapse = " ")))
lines <- c(lines, .render(rows, c("grid", "effect", "points", "sizes", "probit_chi2df",
                                  "logit_chi2df", "simulated_rates")), "")

lines <- c(lines, "Summary by grid shape, over the twelve combinations of effect size and target.", "")
rows <- lapply(names(GRIDS), function(g) {
  k <- grid_of == g
  cv <- vapply(curves, function(r) r$grid == g, logical(1))
  c(g, as.character(sum(k)),
    .pct(mean(num("p_err")[k])), .pct(mean(num("l_err")[k])),
    .num(mean(vapply(curves[cv], function(r) r$pchi, numeric(1))), 2),
    .num(mean(vapply(curves[cv], function(r) r$lchi, numeric(1))), 2),
    paste0(sum(flag("p_cov")[k]), "/", sum(k)),
    paste0(sum(flag("l_cov")[k]), "/", sum(k)))
})
lines <- c(lines, .render(rows, c("grid", "checks", "probit_mae%", "logit_mae%",
                                  "probit_chi2df", "logit_chi2df", "probit_cover",
                                  "logit_cover")), "")

nk <- length(checks)
fw <- num("f_width"); pw <- num("p_width")
both <- is.finite(fw)
lines <- c(lines,
  "Over all checks.",
  "",
  sprintf("  checks                              %d", nk),
  sprintf("  probit mean absolute error          %.2f%%", mean(num("p_err"))),
  sprintf("  probit worst absolute error         %.2f%%", max(num("p_err"))),
  sprintf("  logit mean absolute error           %.2f%%", mean(num("l_err"))),
  sprintf("  logit worst absolute error          %.2f%%", max(num("l_err"))),
  sprintf("  delta-method coverage               %d/%d", sum(flag("p_cov")), nk),
  sprintf("  logit delta-method coverage         %d/%d", sum(flag("l_cov")), nk),
  sprintf("  Fieller coverage                    %d/%d", sum(flag("f_cov")), nk),
  sprintf("  Fieller bounded                     %d/%d", sum(both), nk),
  sprintf("  delta-method mean width             %.2f subjects", mean(pw)),
  sprintf("  Fieller mean width                  %.2f subjects", mean(fw[both])),
  sprintf("  delta-method mean width, bounded    %.2f subjects", mean(pw[both])),
  sprintf("  logit solved outside the range      %d/%d", sum(flag("l_out")), nk),
  sprintf("  solves solve_curve() refused        %d", length(refusals)),
  "")
if (length(refusals)) lines <- c(lines, "Refused, and left out of everything above.", "", refusals, "")

con <- file(out_path, open = "wb")           # binary: LF endings on every platform
writeLines(lines, con, sep = "\n")
close(con)

cat("\n", paste(tail(lines, 20), collapse = "\n"), "\n", sep = "")
cat(sprintf("written to %s in %.0f seconds\n", out_path, elapsed))
