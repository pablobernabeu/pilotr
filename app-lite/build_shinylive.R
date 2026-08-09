# Build the serverless (shinylive / webR) static site for the lite app.
# Stages the lite app.R next to the engine + spec-builder sources (so webR can source them),
# then exports a static site to build/shinylive-demo/. Deploy that folder to GitHub Pages.

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
pkg_R <- file.path(here, "..", "r", "pilotr", "R")
stage <- file.path(here, "..", "build", "shinylive-app")
out   <- file.path(here, "..", "build", "shinylive-demo")

if (!requireNamespace("shinylive", quietly = TRUE))
  install.packages("shinylive", repos = "https://cloud.r-project.org")

unlink(stage, recursive = TRUE); dir.create(stage, recursive = TRUE, showWarnings = FALSE)
# engine-files.txt is the single place the browser build's engine subset is named;
# app.R loads whatever ends up staged, so nothing has to be kept in step by hand.
engine <- readLines(file.path(here, "engine-files.txt"))
engine <- trimws(sub("#.*$", "", engine))
engine <- engine[nzchar(engine)]
stopifnot(all(file.exists(file.path(pkg_R, engine))))
file.copy(file.path(pkg_R, engine), stage)
file.copy(file.path(here, "app.R"), stage, overwrite = TRUE)

cat("Exporting shinylive site (downloads webR assets on first run)...\n")
shinylive::export(stage, out)
cat("Done. Static site in:", normalizePath(out), "\n")
cat("Preview locally with:  httpuv::runStaticServer('", normalizePath(out), "')\n", sep = "")
