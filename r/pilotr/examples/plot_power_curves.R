# Build the publication-style figure: two simulation-based power curves plus the Type M
# (exaggeration) story. Reads the CSVs written by power_curves.py and power_curve_mixed.R.

library(ggplot2)

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
BUILD <- file.path(here, "..", "..", "..", "build")

g <- read.csv(file.path(BUILD, "power_curve_gaussian.csv"))
m_path <- file.path(BUILD, "power_curve_mixed.csv")
m <- if (file.exists(m_path)) read.csv(m_path) else NULL

blue <- "#2C6FB0"; red <- "#B0402C"; grey <- "#888888"
theme_set(theme_minimal(base_size = 12))

# Panel 1: Gaussian power vs N
p1 <- ggplot(g, aes(n_subject, power)) +
  geom_hline(yintercept = 0.8, linetype = 2, colour = grey) +
  geom_line(colour = blue, linewidth = 0.8) + geom_point(colour = blue, size = 2) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = expression("Total " * italic(N)), y = "Power", title = "Two-group Gaussian (d = 0.5)")

# Panel 2: crossed mixed-effects power vs subjects
p2 <- if (!is.null(m)) {
  ggplot(m, aes(n_subject, power)) +
    geom_hline(yintercept = 0.8, linetype = 2, colour = grey) +
    geom_line(colour = red, linewidth = 0.8) + geom_point(colour = red, size = 2) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = expression(italic(N) ~ "subjects (24 items)"), y = "Power",
         title = "Crossed mixed-effects RT")
} else {
  ggplot() +
    annotate("text", x = 0, y = 0, label = "Crossed mixed-effects RT\n(run power_curve_mixed.R)") +
    theme_void()
}

# Panel 3: Type M exaggeration vs power (Gaussian)
p3 <- ggplot(g, aes(power, type_m)) +
  geom_hline(yintercept = 1, linetype = 2, colour = grey) +
  geom_line(colour = blue, linewidth = 0.8) + geom_point(colour = blue, size = 2) +
  labs(x = "Power", y = "Type M (exaggeration ratio)", title = "Low power inflates effects")

fig <- if (requireNamespace("patchwork", quietly = TRUE)) {
  patchwork::wrap_plots(p1, p2, p3, nrow = 1)
} else {
  message("Install 'patchwork' for the combined 3-panel layout; saving panel 1 only.")
  p1
}

ggsave(file.path(BUILD, "power_curves.png"), fig, width = 11.5, height = 3.6, dpi = 130)
cat("wrote", normalizePath(file.path(BUILD, "power_curves.png")), "\n")
