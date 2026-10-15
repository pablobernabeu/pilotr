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
# Without template_params the page is titled "Shiny App" and has no icon, which is what a browser
# tab showed. The site workflow (.github/workflows/site.yml) passes the same values, and the
# icons it links are the R site's, copied beside the page.
icons <- paste(
  '<link rel="icon" type="image/png" sizes="96x96" href="favicon-96x96.png">',
  '<link rel="icon" type="image/svg+xml" href="favicon.svg">', sep = "\n")
shinylive::export(stage, out,
                  template_params = list(title = "pilotr: no-code app", include_in_head = icons))
favicons <- file.path(here, "..", "r", "pilotr", "pkgdown", "favicon")
stopifnot(all(file.copy(file.path(favicons, c("favicon-96x96.png", "favicon.svg")), out,
                        overwrite = TRUE)))
cat("Done. Static site in:", normalizePath(out), "\n")
cat("Preview locally with:  httpuv::runStaticServer('", normalizePath(out), "')\n", sep = "")
