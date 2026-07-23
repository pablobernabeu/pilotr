# Package-level documentation. The "_PACKAGE" sentinel tells roxygen2 to build the
# ?pilotr landing topic, taking the title, description, author and useful links from
# DESCRIPTION, so only the orientation below is written by hand. The topic is
# deliberately not marked @keywords internal: pkgdown's default reference index drops
# internal topics, and this overview is meant to head that index rather than hide from it.

#' @description
#' Generative simulation of experimental and behavioural data sets from a
#' portable JSON design specification shared with the Python package of the
#' same name. Fixed effect sizes are user-specified, by-subject and by-item
#' random intercepts and slopes are crossed, and the response families cover
#' Gaussian, lognormal, shifted lognormal, Bernoulli, Poisson, ordinal and
#' Beta outcomes. Power and precision-based design analysis run from the same
#' specification.
#'
#' @details
#' A pilotr workflow begins with a design specification, a plain list recording
#' the study you plan to run: its groups and conditions, sample sizes, fixed
#' effect sizes, random-effect standard deviations, and the response family.
#' Assemble one from a flat list of design inputs with
#' [build_spec()], or read one back from a JSON file with
#' [load_spec()]. The package ships one ready-to-run specification per
#' design family, and [pilotr_example()] returns their paths.
#' [default_response_name()] gives the response column that a
#' family uses by default, and [spec_json()] serialises a specification
#' back to JSON for the Python twin or the no-code app to read.
#'
#' [simulate_design()] turns a specification into an analysis-ready
#' data frame with one row per observation. A specification carries its own
#' seed, and both languages draw from the combined generator built by
#' [make_rng()] on top of the inverse-normal routine
#' [as241()], so a given specification and seed produce identical data
#' in either language.
#'
#' For analysis, [model_data()] adds the response column that the
#' model expects, and [model_formula()] derives the maximal
#' mixed-model formula the design implies. [brms_bridge()] returns
#' the formula, family and priors for a Bayesian fit.
#'
#' Design analysis runs from that same specification.
#' [power_design()] estimates power for a two-group Gaussian
#' design, together with the Type S and Type M errors of Gelman and Carlin
#' (2014). [power_mixed()] does the same for a crossed
#' mixed-effects design, and [power_curve_mixed()] sweeps
#' sample size to locate where a design becomes adequately powered.
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
#' @author Pablo Bernabeu, author and maintainer
#'   (\email{pcbernabeu@gmail.com},
#'   \href{https://orcid.org/0000-0003-1083-2460}{ORCID}).
"_PACKAGE"
