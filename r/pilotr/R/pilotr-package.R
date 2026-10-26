# Package-level documentation. The "_PACKAGE" sentinel tells roxygen2 to build the
# ?pilotr landing topic, taking the title, description, author and useful links from
# DESCRIPTION, so only the orientation below is written by hand. The topic is
# deliberately not marked @keywords internal: pkgdown's default reference index drops
# internal topics, and this overview is meant to appear in the index, as its last entry.

#' @description
#' Generative simulation of experimental and behavioural data sets from a
#' portable JSON design specification shared with the Python package of the
#' same name. Fixed effect sizes are user-specified, by-subject and by-item
#' random intercepts and slopes are crossed, and the response families cover
#' Gaussian, lognormal, shifted lognormal, ex-Gaussian, Bernoulli, Poisson,
#' ordinal and Beta outcomes. Power and precision-based design analysis run
#' from the same specification.
#'
#' @details
#' A pilotr workflow begins with a design specification, a plain list recording
#' the study you plan to run: its groups and conditions, sample sizes, fixed
#' effect sizes, random-effect standard deviations, and the response family.
#' Assemble one from a flat list of design inputs with
#' [build_spec()], or read one back from a JSON file with
#' [load_spec()]. The package ships eight ready-to-run specifications, and
#' [pilotr_example()] returns their paths.
#' [default_response_name()] gives the response column that a
#' family uses by default, and [spec_json()] serialises a specification
#' back to JSON for the Python twin or the no-code app to read.
#'
#' [simulate_design()] turns a specification into an analysis-ready
#' data frame with one row per observation. A specification carries its own
#' seed, and both languages draw from the combined generator built by
#' [make_rng()] on top of the inverse-normal routine
#' [as241()], so a given specification and seed produce identical data
#' in either language wherever the two share a maths library, as on Linux,
#' outside a rare Beta case under R 4.6.0 or later.
#' Elsewhere, a value that passes through `exp()`, `log()` or `pow()` can
#' differ in its last bits, as `vignette("cross-language", package = "pilotr")`
#' sets out.
#'
#' For analysis, [model_data()] adds the response, contrast and interaction columns that the
#' model expects, and [model_formula()] derives the mixed-model formula
#' the design implies, maximal by default. [brms_bridge()] returns
#' the formula, family and priors for a Bayesian fit.
#'
#' Design analysis runs from that same specification.
#' [power_design()] estimates power for a two-group Gaussian
#' design with one row per subject, together with the Type S and Type M errors
#' of Gelman and Carlin (2014). [power_mixed()] does the same for a crossed
#' mixed-effects design, and [power_curve_mixed()] sweeps
#' sample size. [target_n()] solves that curve for the sample size at which
#' a design becomes adequately powered, with an interval on it.
#' [precision_design()] and its curve counterpart
#' [precision_curve()] report the width of the interval a design
#' buys and the decision probabilities against a region of practical
#' equivalence.
#'
#' Two functions round the package off. [generate_r_script()]
#' writes a self-contained script that reproduces a simulation, and
#' [run_app()] launches the bundled no-code app.
#'
#' For a worked introduction, see
#' `vignette("getting-started", package = "pilotr")`.
#'
#' @references Gelman, A. and Carlin, J. (2014). Beyond power calculations: Assessing Type S
#'   (sign) and Type M (magnitude) errors. \emph{Perspectives on Psychological Science},
#'   9(6), 641-651. \doi{10.1177/1745691614551642}
#'
#' @author Pablo Bernabeu, author and maintainer
#'   (\email{pcbernabeu@gmail.com},
#'   \href{https://orcid.org/0000-0003-1083-2460}{ORCID}).
"_PACKAGE"
