# Build the R site's favicons into r/pilotr/pkgdown/favicon/: the files pkgdown 2.1.2 and later
# link from every page (favicon.ico, favicon-96x96.png, favicon.svg, apple-touch-icon.png and
# site.webmanifest) and the two icons the manifest names. pkgdown copies the folder into the site
# as it stands, so a file missing from it leaves its link returning 404 on every page.
#
# pkgdown::build_favicons() makes the same set by sending the package logo to
# realfavicongenerator.net. This script makes it offline, and from a different master:
# python/docs/assets/favicon.svg, the logo's hexagon and power curve without the wordmark. The
# wordmark cannot be read at 16 or 32 pixels, which is why the R site's favicons, the
# Python site's and both apps' show this design, while build_favicons() would take
# man/figures/logo.svg, wordmark included.
#
# Run from anywhere:  Rscript tools/favicons/build_favicons.R
# It needs the rsvg, magick and jsonlite packages, and it replaces the folder's contents.

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
root <- normalizePath(file.path(here, "..", ".."))
master <- file.path(root, "python", "docs", "assets", "favicon.svg")
out <- file.path(root, "r", "pilotr", "pkgdown", "favicon")

for (pkg in c("rsvg", "magick", "jsonlite"))
  if (!requireNamespace(pkg, quietly = TRUE)) stop("This script needs the ", pkg, " package.")
stopifnot(file.exists(master))

# The folder is emptied first. Writing over it would leave behind any file of an older set that
# the current one does not include.
unlink(out, recursive = TRUE)
dir.create(out, recursive = TRUE)

# The master is square, so each icon is the whole drawing, scaled, on a transparent ground.
png_icon <- function(size, file)
  rsvg::rsvg_png(master, file.path(out, file), width = size, height = size)
png_icon(96, "favicon-96x96.png")
png_icon(180, "apple-touch-icon.png")
png_icon(192, "web-app-manifest-192x192.png")
png_icon(512, "web-app-manifest-512x512.png")

# One .ico holding the three sizes a browser picks from.
ico <- magick::image_join(lapply(c(16, 32, 48), function(size)
  magick::image_read(rsvg::rsvg_png(master, width = size, height = size))))
magick::image_write(ico, file.path(out, "favicon.ico"), format = "ico")

stopifnot(file.copy(master, file.path(out, "favicon.svg")))

# The icon paths are relative, so that they resolve beside the manifest under /pilotr/r/. The
# absolute paths realfavicongenerator.net writes would point at the root of the domain.
manifest <- list(
  name = "pilotr",
  short_name = "pilotr",
  icons = list(
    list(src = "web-app-manifest-192x192.png", sizes = "192x192", type = "image/png"),
    list(src = "web-app-manifest-512x512.png", sizes = "512x512", type = "image/png")),
  theme_color = "#ffffff",
  background_color = "#ffffff",
  display = "standalone")
writeLines(jsonlite::toJSON(manifest, auto_unbox = TRUE, pretty = TRUE),
           file.path(out, "site.webmanifest"), useBytes = TRUE)

cat("Wrote", paste(sort(list.files(out)), collapse = ", "), "to", out, "\n")
