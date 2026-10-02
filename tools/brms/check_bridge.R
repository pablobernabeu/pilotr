# Check that brms accepts the models brms_bridge() and generate_design_analysis() emit, for
# every shipped example and for the two response families no example uses.
#
# Both functions write brms code without loading brms, which the package keeps out of its
# dependencies, Suggests included (design_analysis.R gives the reasons). The package's own tests
# therefore cannot ask brms anything, and this script does. For each design, it simulates a
# data set, adds the analysis columns with model_data() as the emitted code tells its reader to,
# and has brms check the priors against the model and write the Stan program. Three things are
# checked:
#
#   1. The bridge's code, evaluated exactly as emitted, with brm() standing in for a function
#      that runs brms::validate_prior() and brms::stancode() on the arguments it receives.
#   2. The design analysis's brm() call, taken from the emitted script, which adds
#      sample_prior = "yes" and so asks brms to sample from every prior as well.
#   3. Every focal effect the script would test with hypothesis(), which has to be a
#      coefficient brms gives the model, under the name the script reads its draws by.
#
# The designs are the shipped examples and crossed_mixed_rt with a bernoulli and with an
# exgaussian response. No example uses either family, and exgaussian is one of the families
# whose priors are scaled to the data.
#
# Nothing is compiled or sampled, so no C++ toolchain is needed. The whole check takes a few
# minutes, most of them on the larger crossed designs.
#
# Usage: Rscript tools/brms/check_bridge.R
# The exit status is 1 when brms refuses anything, and the output names the design and why.

root <- normalizePath(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "..", ".."), mustWork = FALSE)
if (!nzchar(root) || is.na(root)) root <- normalizePath(".")

if (!requireNamespace("brms", quietly = TRUE))
  stop("this check needs brms, which is not installed", call. = FALSE)

# Source the whole package, as tools/parity/run_r.R does, so that the check runs against the
# code in this checkout.
src <- file.path(root, "r", "pilotr", "R")
for (f in sort(list.files(src, pattern = "\\.R$", full.names = TRUE))) source(f)

# Attached, as the emitted code attaches it, so that prior() and the family constructors in that
# code resolve to brms's own.
suppressPackageStartupMessages(library(brms))
cat("brms", as.character(utils::packageVersion("brms")), "\n")

# The examples pilotr_example() lists, read from the source tree, since system.file() finds
# nothing for a package that is only sourced.
examples <- sort(list.files(file.path(root, "r", "pilotr", "inst", "examples"),
                            pattern = "\\.json$", full.names = TRUE))
if (!length(examples)) stop("no example specifications found", call. = FALSE)

# The crossed example with its response swapped for another family, and its intercept moved to a
# value on that family's scale. It declares specification version 0.3, which the exgaussian
# family needs, and keeps the correlated slopes and so the LKJ prior.
swap_response <- function(response, intercept) {
  spec <- load_spec(file.path(root, "r", "pilotr", "inst", "examples", "crossed_mixed_rt.json"))
  spec[["spec_version"]] <- "0.3"
  spec[["response"]] <- response
  spec[["fixed"]][["intercept"]] <- intercept
  validate_spec(spec)
  spec
}

# Each design is read inside the loop below, so that one that cannot be read is reported as a
# failure of its own.
designs <- c(
  stats::setNames(lapply(examples, function(path) function() load_spec(path)),
                  sub("\\.json$", "", basename(examples))),
  list(crossed_mixed_rt_bernoulli = function()
         swap_response(list(family = "bernoulli", name = "correct"), 1),
       crossed_mixed_rt_exgaussian = function()
         swap_response(list(family = "exgaussian", name = "RT", sigma = 40, beta = 80), 600)))

# Stands in for brms::brm(), taking the arguments the emitted calls pass. It checks the priors
# against the model and writes the Stan program, and fits nothing. The sampler settings arrive
# in `...` and are ignored.
checking_brm <- function(formula, data, family, prior = NULL, sample_prior = "no", ...) {
  brms::validate_prior(prior, formula, data = data, family = family,
                       sample_prior = sample_prior)
  brms::stancode(formula, data = data, family = family, prior = prior,
                 sample_prior = sample_prior)
  invisible(NULL)
}

# The value of the top-level assignment to `name` in a parsed script.
assigned <- function(exprs, name) {
  hit <- Filter(function(e) is.call(e) && identical(e[[1]], as.name("<-")) &&
                  identical(e[[2]], as.name(name)), as.list(exprs))
  if (length(hit) != 1L) stop("the emitted script has no single assignment to `", name, "`")
  hit[[1]][[3]]
}

check_design <- function(spec) {
  bridge <- brms_bridge(spec)
  data <- model_data(spec, simulate_design(spec))

  # 1. The bridge's code as emitted, whose `data` argument is `your_data`.
  env <- new.env()
  env$brm <- checking_brm
  env$your_data <- data
  eval(parse(text = bridge$code), envir = env)

  # 2. The design analysis's call, with the variables the script defines before it.
  keys <- names(spec$fixed$coefficients)
  if (!length(keys)) return(invisible(NULL))
  exprs <- parse(text = generate_design_analysis(spec, focal = keys), keep.source = FALSE)
  env <- new.env()
  env$brm <- checking_brm
  env$dat <- data
  env$spec <- spec
  eval(assigned(exprs, "fit"), envir = env)

  # 3. The script reads each focal effect's draws as b_<name>, so every name has to be one of
  # the model's population-level coefficients.
  focal <- eval(assigned(exprs, "focal"))
  bp <- brms::default_prior(stats::as.formula(bridge$formula), data = data,
                            family = eval(str2lang(bridge$family)))
  coefs <- bp$coef[bp$class == "b" & nzchar(bp$coef)]
  missing <- setdiff(focal, coefs)
  if (length(missing))
    stop("the design analysis would test ", paste(sprintf("'%s'", missing), collapse = ", "),
         ", but brms names the model's coefficients ",
         paste(sprintf("'%s'", coefs), collapse = ", "))
  invisible(NULL)
}

failed <- character(0)
for (name in names(designs)) {
  problem <- tryCatch({
    check_design(designs[[name]]())
    NULL
  }, error = function(e) conditionMessage(e))
  if (is.null(problem)) {
    cat(sprintf("ok    %s\n", name))
  } else {
    cat(sprintf("FAIL  %s: %s\n", name, gsub("\\s+", " ", problem)))
    failed <- c(failed, name)
  }
}

if (length(failed)) {
  cat(sprintf("\nbrms refused %d of %d designs: %s\n", length(failed), length(designs),
              paste(failed, collapse = ", ")))
  quit(status = 1L)
}
cat(sprintf("\nbrms accepted all %d designs\n", length(designs)))
