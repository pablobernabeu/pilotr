# Build the publication-style figure: two simulation-based power curves plus the Type M
# (exaggeration) story. Reads the CSVs written by power_curves.py and power_curve_mixed.R.

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
BUILD <- file.path(here, "..", "..", "..", "build")

g <- read.csv(file.path(BUILD, "power_curve_gaussian.csv"))
m_path <- file.path(BUILD, "power_curve_mixed.csv")
m <- if (file.exists(m_path)) read.csv(m_path) else NULL

blue <- "#2C6FB0"; red <- "#B0402C"; grey <- "#888888"
png(file.path(BUILD, "power_curves.png"), width = 1500, height = 460, res = 130)
par(mfrow = c(1, 3), mar = c(4.2, 4.2, 3, 1), cex.lab = 1.05)

# Panel 1: Gaussian power vs N
plot(g$n_subject, g$power, type = "b", pch = 19, col = blue, lwd = 2,
     ylim = c(0, 1), xlab = "Total N", ylab = "Power",
     main = "Two-group Gaussian (d = 0.5)")
abline(h = 0.8, lty = 2, col = grey); text(min(g$n_subject), 0.83, "80%", col = grey, adj = 0)

# Panel 2: crossed mixed-effects power vs subjects
if (!is.null(m)) {
  plot(m$n_subject, m$power, type = "b", pch = 19, col = red, lwd = 2,
       ylim = c(0, 1), xlab = "Number of subjects (24 items)", ylab = "Power",
       main = "Crossed mixed-effects RT")
  abline(h = 0.8, lty = 2, col = grey); text(min(m$n_subject), 0.83, "80%", col = grey, adj = 0)
} else {
  plot.new(); title(main = "Crossed mixed-effects RT\n(run power_curve_mixed.R)")
}

# Panel 3: Type M exaggeration vs power (Gaussian)
plot(g$power, g$type_m, type = "b", pch = 19, col = blue, lwd = 2,
     xlab = "Power", ylab = "Type M (exaggeration ratio)",
     main = "Low power inflates effects")
abline(h = 1, lty = 2, col = grey); text(max(g$power), 1.05, "no bias", col = grey, adj = 1)

dev.off()
cat("wrote", normalizePath(file.path(BUILD, "power_curves.png")), "\n")
