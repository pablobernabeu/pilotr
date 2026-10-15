# pilotr no-code app. This is the third interface over the shared design spec.
#
# A thin client. Every control writes into the portable JSON spec, which downloads and runs
# unchanged in the R and Python packages. Launch with pilotr::run_app() (installed), or with
# shiny::runApp("r/pilotr/inst/app") (from source).

library(shiny)
library(ggplot2)

# shiny sources this file into an environment whose parent is the global environment, so the
# package's functions are visible here only when pilotr is attached. pilotr::run_app() loads
# the namespace without attaching it, and shiny::runApp() on the installed app directory does
# not attach pilotr either. When this file is the installed package's own app.R, the app
# therefore binds every export of the namespace into its environment, which leaves the search
# path as it was. Anything the app uses from pilotr must accordingly be exported, and helpers
# that only the app needs are defined below. From the source tree, the app sources the engine
# beside it instead, even in a session where an installed pilotr is loaded.
ENGINE_FILES <- NULL   # resolved source paths (dev), which Verify's clean R session sources too
if (!exists("simulate_design", mode = "function")) {
  .installed_app <- system.file("app", package = "pilotr")
  if (nzchar(.installed_app) &&
      identical(normalizePath(getwd(), winslash = "/"),
                normalizePath(.installed_app, winslash = "/")) &&
      requireNamespace("pilotr", quietly = TRUE)) {
    for (.name in getNamespaceExports("pilotr"))
      assign(.name, getExportedValue("pilotr", .name), envir = environment())
  } else {
    # Load the whole engine directory rather than naming its files. The list this
    # replaced had to track R/ by hand, and when validate.R arrived with 0.3.0 it
    # did not: .SPEC_VERSION went undefined and every Simulate click failed. Each
    # file in R/ only defines functions and constants, so order does not matter.
    # "." last, for a layout that stages the engine beside the app; app.R and the
    # roxygen sentinel are dropped so that case cannot source this file into itself.
    # A directory counts only if it holds simulate.R. An installed package has an R/
    # directory too, holding the lazy-load database and no sources, and accepting any
    # directory sent the installed app there.
    .dir <- Find(function(d) file.exists(file.path(d, "simulate.R")),
                 c("../../R", "../R", "R", "."))
    if (is.null(.dir))
      stop("cannot find the pilotr engine. Launch the app with pilotr::run_app(), or run it ",
           "from the package's source tree, whose R/ directory holds simulate.R.")
    .files <- sort(list.files(.dir, pattern = "[.][Rr]$", full.names = TRUE))
    .files <- .files[!basename(.files) %in% c("pilotr-package.R", "app.R")]
    for (f in .files) source(f)
    ENGINE_FILES <- normalizePath(.files)
  }
}

MAX_SIMS <- as.integer(Sys.getenv("PILOTR_MAX_SIMS", "5000"))
N_SIMS_DEFAULT <- 1000L   # the value the Simulations box starts at
N_SIMS_MIN     <- 100L    # the smallest count the box accepts

# The replicate count to run, read from the Simulations box and clamped to the app's range.
# A cleared box arrives as NA and a pasted word arrives as NA with a warning, and clamping NA
# leaves NA: min(max(NA, 100), 5000) is NA, and power_design() then stops with "vector size
# cannot be NA", which tells the user nothing about the box they emptied. The fallback is the
# value the box started at, which is the count they would have run anyway.
.n_sims_input <- function(v) {
  v <- suppressWarnings(as.integer(v))
  if (length(v) != 1L || is.na(v)) N_SIMS_DEFAULT else min(max(v, N_SIMS_MIN), MAX_SIMS)
}

# Type S and Type M are shown with the Monte Carlo error power_design() reports for each. Type S
# has a Wilson interval and Type M a standard error. A rate whose error is undefined is shown
# bare, as before. That happens to both at a zero true effect, and to Type M when a single
# replicate is significant. app-lite/app.R formats them in the same way.
type_s_text <- function(r)
  if (is.na(r$type_s_lo)) sprintf("%.4f", r$type_s) else
    sprintf("%.4f (95%% Wilson [%.4f, %.4f])", r$type_s, r$type_s_lo, r$type_s_hi)
type_m_text <- function(r)
  if (is.na(r$type_m_mcse)) sprintf("%.3f", r$type_m) else
    sprintf("%.3f ± %.3f (MCSE)", r$type_m, r$type_m_mcse)

# Where the data Verify's script produced first differ from the reference data, as the phrase
# its verdict ends with, or NULL when the two are identical. Doubles are compared bit by bit
# (num.eq = FALSE), so that a zero of the other sign counts as a difference, as the verdict's
# "bit for bit" says. Verify used to compare the sums of the response column with a tolerance,
# which let through reversed rows and any change to the other columns.
.data_difference <- function(ref, script_data) {
  if (!identical(names(script_data), names(ref)))
    return(sprintf("the column names: the script gives %s, and the specification %s",
                   paste(names(script_data), collapse = ", "), paste(names(ref), collapse = ", ")))
  if (nrow(script_data) != nrow(ref))
    return(sprintf("the number of rows: the script gives %d, and the specification %d",
                   nrow(script_data), nrow(ref)))
  for (column in names(ref)) {
    a <- ref[[column]]
    b <- script_data[[column]]
    if (identical(a, b, num.eq = FALSE)) next
    if (!identical(typeof(a), typeof(b)) || !identical(attributes(a), attributes(b)))
      return(sprintf("the type or attributes of column '%s'", column))
    a <- unclass(a)
    b <- unclass(b)
    first_row <- Find(function(i) !identical(a[[i]], b[[i]], num.eq = FALSE), seq_along(a))
    return(if (is.null(first_row)) sprintf("column '%s'", column)
           else sprintf("column '%s', at row %d", column, first_row))
  }
  if (identical(ref, script_data, num.eq = FALSE)) NULL
  else "the attributes of the data frame, such as its row names"
}

# Power runs in a background R process only when run_app() asks for it, through this option, and
# only from the installed package, which a worker can load. From source, and under
# shiny::runApp() on this directory, power runs in this process. The app used to run power in
# the background whenever future and promises were installed, so run_app(async = FALSE) changed
# nothing. A background run leaves this process free for other sessions. The session that started
# it still waits for the result, since Shiny holds a session's inputs and outputs until the
# promise its observer returned has settled.
.async_ok <- isTRUE(getOption("pilotr.app_async")) && isNamespaceLoaded("pilotr") &&
  nzchar(system.file(package = "future")) && nzchar(system.file(package = "promises"))

# The power curve over N subjects, on a grid extended until it brackets the 0.80 target. The grid
# starts at half, three quarters, one, one and a half and twice `n`, at least 4. While the top
# point's power falls short of the target by more than twice its Monte Carlo standard error,
# sqrt(p (1 - p) / n_sims), a point half as large again is added above it, up to `cap` subjects.
# While the bottom point's power exceeds the target by more than that, a point two thirds its
# size is added below it, down to 4. The grid holds one point for each of `seeds`, so 10 at most.
# The fixed grid stopped at twice N, which for the default design is the target itself (127.5
# subjects by power.t.test()). About half of all runs were then refused, and the rest solved low,
# since a run solved only when its top point happened to land above 0.80.
#
# power_at(size, seed) gives the power at `size` subjects. Each point takes the next of `seeds` in
# the order the points are computed, so a point added below the grid never changes the seed, or
# the power, of a point already computed. Returns the grid in increasing order, with the power and
# the seed at each point. app-lite/app.R carries the same function, with a lower cap.
power_curve_grid <- function(n, cap, n_sims, seeds, power_at, target = 0.8) {
  grid <- unique(pmax(4, round(n * c(0.5, 0.75, 1, 1.5, 2))))
  pw <- vapply(seq_along(grid), function(k) power_at(grid[k], seeds[k]), numeric(1))
  margin <- function(p) 2 * sqrt(p * (1 - p) / n_sims)
  while (length(grid) < length(seeds)) {
    top <- max(grid)
    p <- pw[which.max(grid)]
    if (top >= cap || p >= target + margin(p)) break
    size <- min(cap, round(top * 1.5))
    pw <- c(pw, power_at(size, seeds[length(grid) + 1L]))
    grid <- c(grid, size)
  }
  while (length(grid) < length(seeds)) {
    bottom <- min(grid)
    p <- pw[which.min(grid)]
    if (bottom <= 4 || p <= target - margin(p)) break
    size <- max(4, round(bottom / 1.5))
    pw <- c(pw, power_at(size, seeds[length(grid) + 1L]))
    grid <- c(grid, size)
  }
  ord <- order(grid)
  list(grid = grid[ord], pw = pw[ord], seeds = seeds[seq_along(grid)][ord])
}

# The total N at which the two-sample t-test that power_design() simulates reaches a power of 0.80,
# by stats::power.t.test(), shown beside the curve's own solve as a check on it. Only for a plain
# design, without random effects or continuous predictors, since either adds variance that the
# `sd` of power.t.test() leaves out. NULL otherwise, and for an effect of zero, which no sample
# size detects. app-lite/app.R carries the same function.
t_test_n80 <- function(spec) {
  if (length(spec[["random"]]) || length(spec[["predictors"]])) return(NULL)
  tryCatch({
    f <- Filter(function(x) !is.null(x[["between"]]), spec[["factors"]])[[1]]
    contrast <- names(f[["contrasts"]])[1]
    delta <- abs(spec[["fixed"]][["coefficients"]][[contrast]] *
                   diff(f[["contrasts"]][[contrast]]))
    if (length(delta) != 1L || !is.finite(delta) || delta == 0) NULL
    else 2 * stats::power.t.test(delta = delta, sd = spec[["response"]][["sigma"]],
                                 power = 0.8)$n
  }, error = function(e) NULL)
}

# Reasonable intercept/effect (and noise) defaults per family, on each family's own scale.
# Switching family resets these so a point-and-click design stays valid: a Gaussian-scale
# intercept of 100, for instance, would overflow a log or logit family.
FAMILY_DEFAULTS <- list(
  gaussian          = list(intercept = 100, effect = 5,   sigma = 10),
  shifted_lognormal = list(intercept = 6,   effect = 0.1, sigma = 0.3),
  bernoulli         = list(intercept = 0,   effect = 0.5),
  poisson           = list(intercept = 1.5, effect = 0.3),
  ordinal           = list(intercept = 0,   effect = 0.8),
  beta              = list(intercept = 0,   effect = 0.8, phi = 8)
)

# Fill colours for the plots, as many as the data has levels, from the palette of app-lite/app.R. A
# pasted specification can carry a factor with any number of levels, and the fixed pair of colours
# used before stopped the plot of a three-level design with "Insufficient values in manual scale".
PALETTE <- c("#2C6FB0", "#B0402C", "#2E8B57", "#8E6FB0", "#C8922A", "#5A5A5A")
pal <- function(k) if (k <= length(PALETTE)) PALETTE[seq_len(k)] else grDevices::hcl.colors(k, "Dynamic")

# Shown above the tabs while the controls differ from the specification of the last simulation, in
# the words app-lite/app.R uses.
CONTROLS_CHANGED <- paste0(
  "The controls have changed since the last Simulate: the JSON and the specification download ",
  "follow the controls, while the data, summary, plot and data download show the last ",
  "simulation. Select Simulate to bring them together.")

guide_tab <- tabPanel(
  "Guide",
  tags$div(
    style = "max-width: 52rem;",
    tags$p(tags$b("pilotr"), " builds a portable design specification and simulates ",
           "experimental and behavioural data from it. This is the no-code interface. The ",
           "same specification runs unchanged in the R and Python packages and reproduces ",
           "the same data, bit for bit wherever the two share a maths library. Elsewhere, a ",
           "value that passes through ", tags$code("exp()"), ", ", tags$code("log()"), " or ",
           tags$code("pow()"), " can differ in its last bits, as the ",
           tags$a(href = "https://pablobernabeu.github.io/pilotr/r/articles/cross-language.html",
                  target = "_blank", "cross-language guide"), " explains."),
    tags$ol(
      tags$li("Describe the design on the left: the sample sizes, the factor and its two ",
              "levels, the fixed intercept and effect, and a response family. Changing the ",
              "family resets the intercept and effect to sensible values for that family's ",
              "scale. For within-subjects and crossed designs, add by-subject and by-item ",
              "random effects. ", tags$b("Counterbalanced lists"), " show each subject each item ",
              "in one condition, rotating the conditions over subjects and items, which halves ",
              "the trials of the crossed design."),
      tags$li("Select ", tags$b("Simulate"), ", then read the tabs: the JSON specification, ",
              "the simulated data, a summary and plot, a simulation-based power analysis, and ",
              "a reproducible R script that you can verify in a clean R session."),
      tags$li("Download the specification (", tags$code(".json"), ") or the data (",
              tags$code(".csv"), ") to take the design into the R or Python package, or onto ",
              "a cluster for large power and precision analyses.")
    ),
    tags$h5("Response families"),
    tags$p("Gaussian, shifted lognormal (reaction times), Bernoulli (accuracy), Poisson ",
           "(counts), ordinal (Likert) and Beta (proportions). The effect is the difference ",
           "between the two levels on the response scale: the identity scale for Gaussian, ",
           "the log scale for reaction times and Poisson, and the logit scale for accuracy, ",
           "ordinal and Beta. The ", tags$b("Advanced: paste a JSON spec"), " box accepts ",
           "designs beyond the point-and-click controls, such as continuous predictors, ",
           "interactions and nesting."),
    tags$h5("Power and design analysis"),
    tags$p("The in-app backend estimates power for the two-group Gaussian design with one row ",
           "per subject. It reports power and the Type S and Type M errors of Gelman and Carlin ",
           "(2014), each with its Monte Carlo error, and draws a power curve over sample size, ",
           "adding sample sizes until the curve brackets a power of 0.80. It then solves that ",
           "curve for the sample size at which power reaches 0.80 and reports a confidence ",
           "interval on it. The heavier analyses stay with the packages themselves: crossed ",
           "mixed-effects power in R (lme4, the reference) and in Python (statsmodels), and ",
           "precision/ROPE design analysis in R. The specification you build here drives all ",
           "three interfaces."),
    tags$p(
      tags$a(href = "https://pablobernabeu.github.io/pilotr/", target = "_blank", "Documentation"), " · ",
      tags$a(href = "https://github.com/pablobernabeu/pilotr", target = "_blank", "Source (R and Python)")
    )
  )
)

# ---------------------------------------------------------------- UI ----
ui <- fluidPage(
  tags$head(tags$link(rel = "icon", type = "image/png", href = "favicon.png")),
  titlePanel("pilotr: design, simulate and power analysis (one spec, three interfaces)"),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      textInput("name", "Design name", "my_design"),
      numericInput("seed", "Seed", 2024, step = 1),
      fluidRow(
        column(6, selectInput("design_kind", "Design",
                              c("Between-subjects" = "between", "Within / crossed" = "within"))),
        column(6, numericInput("n_subject", "N subjects", 64, min = 2))
      ),
      conditionalPanel("input.design_kind == 'within'",
        checkboxInput("include_items", "Crossed with items (by-item random effects)", TRUE),
        conditionalPanel("input.include_items",
          numericInput("n_item", "N items", 24, min = 2),
          checkboxInput("counterbalance",
                        "Counterbalanced lists (each subject sees each item in one condition)",
                        FALSE))
      ),
      tags$hr(),
      fluidRow(
        column(6, textInput("factor_name", "Factor name", "group")),
        column(3, textInput("lev1", "Level 1", "control")),
        column(3, textInput("lev2", "Level 2", "treatment"))
      ),
      fluidRow(
        column(6, numericInput("intercept", "Intercept (mean / log-rate / logit)", 100)),
        column(6, numericInput("effect", "Effect (coef. on -0.5/+0.5 contrast)", 5))
      ),
      conditionalPanel("input.design_kind == 'within'",
        tags$b("By-subject random effects"),
        fluidRow(
          column(4, numericInput("subj_int_sd", "Intercept SD", 0.12, min = 0)),
          column(4, numericInput("subj_slope_sd", "Slope SD", 0.04, min = 0)),
          column(4, numericInput("subj_corr", "corr", 0.2, min = -1, max = 1))
        ),
        conditionalPanel("input.include_items",
          tags$b("By-item random effects"),
          fluidRow(
            column(4, numericInput("item_int_sd", "Intercept SD", 0.08, min = 0)),
            column(4, numericInput("item_slope_sd", "Slope SD", 0.02, min = 0)),
            column(4, numericInput("item_corr", "corr", -0.1, min = -1, max = 1))
          )
        )
      ),
      tags$hr(),
      selectInput("family", "Response family",
                  c("Gaussian" = "gaussian", "Shifted lognormal (RT)" = "shifted_lognormal",
                    "Bernoulli (accuracy, logit)" = "bernoulli", "Poisson (counts, log)" = "poisson",
                    "Ordinal (Likert, cumulative-logit)" = "ordinal",
                    "Beta (proportion, logit)" = "beta")),
      textInput("resp_name", "Response name (blank = auto)", ""),
      conditionalPanel("input.family == 'gaussian' || input.family == 'shifted_lognormal'",
        numericInput("sigma", "Residual SD (log scale for RT)", 10, min = 0)),
      conditionalPanel("input.family == 'shifted_lognormal'",
        numericInput("shift", "Shift / non-decision time", 200)),
      conditionalPanel("input.family == 'ordinal'",
        textInput("thresholds", "Thresholds (comma-separated)", "-2, -0.6, 0.6, 2")),
      conditionalPanel("input.family == 'beta'",
        numericInput("phi", "Precision (phi)", 8, min = 0.1)),
      tags$hr(),
      tags$details(
        tags$summary("Advanced: paste a JSON spec (overrides the controls)"),
        # A label, as in the browser build. A NULL one renders an empty <label>, which names the
        # box for nobody, and the summary above is not its label, so a screen reader had only
        # the placeholder to read out.
        textAreaInput("spec_json_in", "Design spec (JSON)", "", rows = 5,
          placeholder = "Paste a pilotr spec with continuous predictors / interactions (e.g. a reading-time design)")),
      actionButton("simulate", "Simulate", class = "btn-primary"),
      downloadButton("dl_spec", "Download spec (.json)"),
      downloadButton("dl_data", "Download data (.csv)")
    ),
    mainPanel(
      width = 8,
      uiOutput("sim_note"),
      tabsetPanel(
        id = "tabs",
        guide_tab,
        tabPanel("Design spec (JSON)",
          p("This portable spec is the single source of truth. Download it and run it ",
            "unchanged in R or Python to obtain the same data set."),
          verbatimTextOutput("json")),
        tabPanel("Data", verbatimTextOutput("dims"), tableOutput("head")),
        tabPanel("Summary & plot", verbatimTextOutput("summary"), plotOutput("plot", height = "320px")),
        tabPanel("Power & design analysis",
          p("Simulation-based power with Type S / Type M (Gelman & Carlin, 2014), for the ",
            "two-group Gaussian design with one row per subject."),
          numericInput("n_sims", "Simulations (capped in-app)", N_SIMS_DEFAULT,
                       min = N_SIMS_MIN, max = MAX_SIMS, step = 100),
          div(class = "mb-3",
              actionButton("run_power", "Run power analysis", class = "btn-primary"), " ",
              actionButton("run_curve", "Power curve over N")),
          verbatimTextOutput("power_out"),
          plotOutput("power_plot", height = "300px")),
        tabPanel("Reproducible R script",
          p("Your no-code design as a self-contained R script. Download it, or select Verify to ",
            "run it in a clean R session and check that its data are identical, bit for bit, to ",
            "the data the current specification produces."),
          downloadButton("dl_rscript", "Download .R"),
          actionButton("verify_code", "Verify (clean R session)"),
          verbatimTextOutput("verify_out"),
          tags$hr(),
          verbatimTextOutput("rscript"),
          tags$hr(),
          verbatimTextOutput("repro_py"))
      )
    )
  )
)

# ------------------------------------------------------------ server ----
server <- function(input, output, session) {

  .grp_col <- function(spec) if (length(spec$factors)) spec$factors[[1]]$name else NULL

  # Reset intercept/effect (and noise) to family-appropriate values when the family changes,
  # so a point-and-click design stays on the right scale (no exp(100) for a log family).
  observeEvent(input$family, {
    d <- FAMILY_DEFAULTS[[input$family]]; if (is.null(d)) return()
    updateNumericInput(session, "intercept", value = d$intercept)
    updateNumericInput(session, "effect",    value = d$effect)
    if (!is.null(d$sigma)) updateNumericInput(session, "sigma", value = d$sigma)
    if (!is.null(d$phi))   updateNumericInput(session, "phi",   value = d$phi)
  }, ignoreInit = TRUE)

  current_spec <- reactive({
    txt <- input$spec_json_in
    if (!is.null(txt) && nzchar(trimws(txt))) {                 # advanced override
      # As in the lite app, route the paste through load_spec() (staged in a temp file) rather
      # than a bare fromJSON, so a pasted spec faces the same validation as a loaded one.
      # Several of the ways a spec can be wrong produce plausible data rather than an error
      # (a mistyped coefficient key resolves to no column and silently generates null-effect
      # data), which is exactly what 0.3.0's validation refuses. validate() surfaces the
      # message in every output that reads the spec.
      spec <- tryCatch({
        tf <- tempfile(fileext = ".json"); writeLines(txt, tf); load_spec(tf)
      }, error = function(e) e)
      if (inherits(spec, "error")) validate(need(FALSE, conditionMessage(spec)))
      return(spec)
    }
    build_spec(list(
      name = input$name, seed = input$seed, n_subject = input$n_subject,
      include_items = input$include_items, n_item = input$n_item,
      counterbalance = input$counterbalance,
      design_kind = input$design_kind, factor_name = input$factor_name,
      lev1 = input$lev1, lev2 = input$lev2, intercept = input$intercept, effect = input$effect,
      subj_int_sd = input$subj_int_sd, subj_slope_sd = input$subj_slope_sd, subj_corr = input$subj_corr,
      item_int_sd = input$item_int_sd, item_slope_sd = input$item_slope_sd, item_corr = input$item_corr,
      family = input$family, resp_name = input$resp_name, sigma = input$sigma,
      shift = input$shift, thresholds = input$thresholds, phi = input$phi
    ))
  })

  # Simulate keeps the specification it ran beside the data it made, and the Data tab, the summary,
  # the plot and the data download read the two together. The summary and the plot used to take
  # the response and factor names from the controls, so a new family or response name after a
  # simulation stopped both. The Data tab printed the Seed box, and the download simulated afresh
  # from the controls. The JSON tab and the specification download follow the controls, and
  # sim_note says when the two have parted.
  sim <- eventReactive(input$simulate, {
    s <- current_spec()
    list(spec = s, data = simulate_design(s))
  }, ignoreNULL = FALSE)

  output$json <- renderText(spec_json(current_spec()))
  # The seed through format(), since %d stops on a whole double beyond R's integer range and
  # build_spec() keeps such a seed.
  output$dims <- renderText({
    s <- sim()
    sprintf("Simulated %d rows x %d columns (seed %s).", nrow(s$data), ncol(s$data),
            format(s$spec$seed, scientific = FALSE))
  })
  output$head <- renderTable(head(sim()$data, 10), striped = TRUE, spacing = "xs")
  # A simulation that failed shows nothing that could be out of date. A specification the controls
  # cannot build counts as a change, since the JSON tab then gives the reason.
  output$sim_note <- renderUI({
    last <- tryCatch(sim(), error = function(e) NULL)
    if (is.null(last)) return(NULL)
    live <- tryCatch(spec_json(current_spec()), error = function(e) NULL)
    if (identical(live, spec_json(last$spec))) return(NULL)
    tags$div(class = "alert alert-warning", role = "status", CONTROLS_CHANGED)
  })

  output$summary <- renderPrint({
    s <- sim(); d <- s$data; spec <- s$spec; yn <- spec$response$name; fn <- .grp_col(spec)
    has_grp <- !is.null(fn) && fn %in% names(d)
    if (is.numeric(d[[yn]])) {
      if (has_grp) {
        agg <- aggregate(d[[yn]], list(d[[fn]]), function(x) c(mean = mean(x), sd = sd(x), n = length(x)))
        cat("Mean (SD) of", yn, "by", fn, ":\n"); print(do.call(data.frame, agg))
      } else cat(sprintf("%s: mean %.4f, SD %.4f, n %d\n", yn, mean(d[[yn]]), sd(d[[yn]]), nrow(d)))
    } else if (has_grp) { cat("Counts of", yn, "by", fn, ":\n"); print(table(d[[fn]], d[[yn]])) }
    else print(table(d[[yn]]))
  })

  output$plot <- renderPlot({
    s <- sim(); d <- s$data; spec <- s$spec; yn <- spec$response$name; fn <- .grp_col(spec)
    has_grp <- !is.null(fn) && fn %in% names(d)
    base <- theme_minimal(base_size = 14)
    if (is.numeric(d[[yn]])) {
      if (has_grp)
        ggplot(d, aes(.data[[fn]], .data[[yn]], fill = .data[[fn]])) +
          geom_boxplot(alpha = 0.85, outlier.alpha = 0.35) +
          scale_fill_manual(values = pal(length(unique(d[[fn]]))), guide = "none") +
          labs(x = fn, y = yn, title = paste("Distribution of", yn)) + base
      else
        ggplot(d, aes(.data[[yn]])) +
          geom_histogram(bins = 30, fill = "#2C6FB0", colour = "white") +
          labs(x = yn, y = "count", title = paste("Distribution of", yn)) + base
    } else if (has_grp)
      ggplot(d, aes(.data[[yn]], fill = .data[[fn]])) +
        geom_bar(position = "dodge") +
        scale_fill_manual(values = pal(length(unique(d[[fn]]))), name = fn) +
        labs(x = yn, y = "count", title = paste("Counts of", yn)) + base
    else
      ggplot(d, aes(.data[[yn]])) +
        geom_bar(fill = "#2C6FB0") +
        labs(x = yn, y = "count", title = paste("Counts of", yn)) + base
  })

  output$rscript <- renderText(generate_r_script(current_spec()))
  output$repro_py <- renderText(paste0(
    "# The same design also runs in Python, with the same data given the same seed:\n",
    "from pilotr import simulate\n",
    "d = simulate(\"design.json\")   # download the spec from the first tab\n",
    "d.to_csv(\"data.csv\")"))

  # Text downloads go through a binary connection so they carry LF endings and a single
  # trailing newline on every platform; text-mode writeLines() writes CRLF on Windows and
  # appends a second newline to text that already ends in one.
  write_text_download <- function(text, file) {
    con <- base::file(file, open = "wb")
    on.exit(close(con), add = TRUE)
    if (!endsWith(text, "\n")) text <- paste0(text, "\n")
    writeLines(text, con, sep = "")
  }

  # A download is named after the specification it holds, with "design" standing in for an empty
  # name, as in the browser build. The specification and script downloads took the Design name
  # box, so a pasted specification was saved under the name in that box, and an emptied box named
  # the files ".json", ".R" and "_seed2024.csv", which are hidden on Unix. A pasted specification
  # that load_spec() refuses has no name to give.
  file_stem <- function(name)
    if (is.character(name) && length(name) == 1L && !is.na(name) && nzchar(name)) name else "design"
  current_stem <- function() file_stem(tryCatch(current_spec()[["name"]], error = function(e) NULL))

  output$dl_rscript <- downloadHandler(
    filename = function() paste0(current_stem(), ".R"),
    content = function(file) write_text_download(generate_r_script(current_spec()), file))

  # ---- Verify: run the exported script in a clean R subprocess and compare ----
  # The clean session runs the text that this tab shows and the download carries, so the verdict
  # covers what the user takes away. It used to simulate from the JSON there and compare the sums
  # of the response column, so it certified a script it never ran. The reference is simulated here
  # from the current specification, which may have changed since the Data tab was last simulated.
  verify_result <- reactiveVal(NULL)
  observeEvent(input$verify_code, {
    spec <- current_spec()
    # A refusal from the engine is reported here. Unguarded, it ended the session, as clearing the
    # N subjects box showed.
    ref <- tryCatch(simulate_design(spec), error = function(e) e)
    if (inherits(ref, "error")) { verify_result(list(msg = paste("error:", conditionMessage(ref)))); return() }
    if (!requireNamespace("callr", quietly = TRUE)) {
      verify_result(list(msg = "Install the 'callr' package to verify in a clean R session.")); return()
    }
    script <- generate_r_script(spec)
    # From the source tree, the script has to run on the engine this app sourced, whatever pilotr
    # may be installed. The clean session sources the same files, and the script's library(pilotr)
    # line is dropped.
    if (!is.null(ENGINE_FILES)) script <- sub("library(pilotr)\n", "", script, fixed = TRUE)
    withProgress(message = "Running the exported script in a clean R session...", value = 0.5, {
      res <- tryCatch(callr::r(function(script, files) {
        for (f in files) source(f)
        env <- new.env(parent = globalenv())
        eval(parse(text = script), envir = env)
        env[["data"]]
      }, args = list(script = script, files = ENGINE_FILES)), error = function(e) e)
    })
    if (inherits(res, "error")) { verify_result(list(msg = paste("error:", conditionMessage(res)))); return() }
    if (!is.data.frame(res)) {
      verify_result(list(msg = "The script ran in a clean R session but left no data frame called data."))
      return()
    }
    verify_result(list(diff = .data_difference(ref, res), n = nrow(ref)))
  })
  output$verify_out <- renderText({
    r <- verify_result()
    if (is.null(r)) return(paste0(
      "Select Verify to run this script in a clean R session and check that its data are ",
      "identical, bit for bit, to the data the current specification produces."))
    if (!is.null(r$msg)) return(r$msg)
    if (is.null(r$diff)) sprintf(paste0(
      "Reproduces bit for bit in a clean R session: the script's data equal the data the ",
      "current specification produces (%d rows, all columns identical)."), r$n)
    else paste0(
      "Does not reproduce in a clean R session: the script's data differ from the data the ",
      "current specification produces. The first difference is in ", r$diff, ".")
  })

  # ---- power: point estimate + curve, capped, in a worker when run_app() asks for it ----
  power_result     <- reactiveVal(NULL)
  power_curve_data <- reactiveVal(NULL)
  # power_design() t-tests every row, so it runs only for a Gaussian design with one two-level
  # factor between subjects, one row per subject and no shared cluster. That rules out an item
  # unit and any grouping factor besides subject, while an item entry groups nothing without an
  # item unit. This is power_design()'s own rule, written out because the installed app reaches
  # only pilotr's exports. Testing the first factor alone let a pasted crossed design through to
  # a t-test of its correlated rows, and let a three-level or 2 x 2 design stop the observer.
  gaussian_two_group <- function(spec) {
    factors <- spec[["factors"]]
    identical(spec[["response"]][["family"]], "gaussian") && length(factors) == 1L &&
      identical(factors[[1]][["between"]], "subject") &&
      is.null(factors[[1]][["vary_within"]]) && length(factors[[1]][["levels"]]) == 2L &&
      is.null(spec[["units"]][["item"]]) &&
      all(names(spec[["random"]]) %in% c("subject", "item"))
  }
  # What to run instead, chosen by what can fit the design. R's power_mixed() fits every family
  # but two on the scale of its coefficients, by lmer() or glmer(), or lm() and glm() without
  # random effects. Python's power_mixed() takes one within factor crossed with items and fits a
  # linear model on the response's own scale. It is offered only for such a design, in a family
  # whose coefficients are on that scale or on the log scale it analyses. Ordinal and Beta
  # outcomes have no frequentist model here and go to the Bayesian design analysis. A single text
  # used to send every design to both power_mixed() functions, including designs neither fits on
  # the scale of their coefficients. The installed app reaches only pilotr's exports, so the rule
  # is written out here, as it is in the browser build.
  not_supported <- function(spec) {
    family <- spec[["response"]][["family"]]
    opening <- paste0(
      "The in-app power backend covers the two-group Gaussian design with one row per\n",
      "subject. ")
    if (family %in% c("ordinal", "beta")) {
      focal <- names(spec[["fixed"]][["coefficients"]])[1]
      if (!length(focal) || is.na(focal)) focal <- "effect"
      return(paste0(opening,
        "pilotr has no frequentist model for ordinal or Beta outcomes yet, so the\n",
        "analysis of this design is Bayesian, on the logit scale its coefficients are\n",
        "written on. Download the spec (the Design spec tab) and run it in R:\n\n",
        "    library(pilotr)\n",
        "    spec <- load_spec(\"design.json\")\n",
        "    brms_bridge(spec)   # the brms model, as code\n",
        sprintf("    generate_design_analysis(spec, focal = \"%s\", file = \"design_analysis.R\")\n",
                focal),
        "\ngenerate_design_analysis() writes a script that simulates the design, fits that\n",
        "model with brms and decides about the effect, replicate by replicate."))
    }
    within <- Filter(function(f) !is.null(f[["vary_within"]]), spec[["factors"]])
    python <- length(within) == 1L && !is.null(spec[["units"]][["item"]]) &&
      family %in% c("gaussian", "lognormal", "shifted_lognormal", "exgaussian")
    paste0(opening,
      "For this design, download the spec (the Design spec tab) and run it\n",
      "directly:\n\n",
      "R (lme4; may take a few minutes):\n",
      "    library(pilotr)\n",
      "    spec <- load_spec(\"design.json\")\n",
      "    power_mixed(spec, n_sims = 200)\n",
      "    power_curve_mixed(spec, subject_ns = c(20, 40, 60), n_sims = 200)",
      if (python) paste0(
        "\n\nPython (statsmodels backend):\n",
        "    from pilotr import load_spec, power_mixed\n",
        "    power_mixed(load_spec(\"design.json\"), n_sims=200)"))
  }

  # Power needs at least 4 subjects, two to a group. Below that, power_design() refuses the design,
  # and the curve keeps too few grid points of 4 or more for target_n() to solve. A cleared N
  # subjects box reaches the specification as NA.
  n_subject_problem <- function(spec) {
    n <- spec[["units"]][["subject"]][["n"]]
    if (length(n) != 1L || is.na(n) || n < 4) "Set N subjects to at least 4 before running power."
  }
  # An error in an observer closes the session, and the design goes with it, so the power buttons
  # report any error in the Power tab. A pasted specification that load_spec() refuses is one. It
  # reaches the observer as a validation error, which Shiny drops without a word when nothing
  # catches it, and the click used to do nothing.
  power_failed <- function(e) power_result(list(msg = paste("error:", conditionMessage(e))))

  observeEvent(input$run_power, {
    power_curve_data(NULL)
    tryCatch({
      spec <- current_spec()
      problem <- if (gaussian_two_group(spec)) n_subject_problem(spec) else not_supported(spec)
      if (!is.null(problem)) { power_result(list(msg = problem)); return() }
      n <- .n_sims_input(input$n_sims)
      if (.async_ok) {
        power_result(list(msg = sprintf("Running %d simulations in a background worker...", n)))
        p <- promises::future_promise({ pilotr::power_design(spec, n_sims = n) }, seed = TRUE)
        promises::then(p, onFulfilled = function(res) power_result(res),
                          onRejected = function(e) power_result(list(msg = paste("error:", conditionMessage(e)))))
      } else {
        withProgress(message = sprintf("Simulating %d datasets...", n), value = 0.5,
                     power_result(power_design(spec, n_sims = n)))
      }
    }, error = power_failed)
  })

  # Reports a computed curve: the text in the Power tab and the curve for the plot.
  report_curve <- function(curve, spec, n) {
    # A dashed target line leaves the reader to judge the crossing. target_n() estimates
    # it instead, and reports its refusal when the curve does not settle the question,
    # which is more use than a number the sweep cannot support.
    solved <- tryCatch(target_n(data.frame(n_subject = curve$grid, power = curve$pw, n_sims = n),
                                target = 0.8),
                       error = function(e) conditionMessage(e))
    reference <- t_test_n80(spec)
    power_result(list(msg = paste0(
      sprintf("Power curve at n_sims = %d per point. N subjects = %s.\n",
              n, paste(curve$grid, collapse = ", ")),
      if (is.character(solved)) paste0("Target power 0.80: ", solved)
      else sprintf("Target power 0.80: N = %d subjects (95%% interval %d to %d).",
                   solved$n, solved$n_lo, solved$n_hi),
      if (!is.null(reference))
        sprintf("\nAnalytic reference: N = %.1f subjects in two equal groups, by stats::power.t.test().",
                reference))))
    power_curve_data(c(curve, list(solved = if (is.character(solved)) NULL else solved)))
  }

  observeEvent(input$run_curve, {
    power_curve_data(NULL)
    tryCatch({
      spec <- current_spec()
      problem <- if (gaussian_two_group(spec)) n_subject_problem(spec) else not_supported(spec)
      if (!is.null(problem)) { power_result(list(msg = problem)); return() }
      n <- .n_sims_input(input$n_sims)
      base_n <- spec$units$subject$n
      cap <- 8 * base_n
      # Each point takes its own seed, so that target_n() fits independent points, as it assumes.
      # Every point used to keep the design's seed.
      seeds <- replicate_seeds(spec$seed, 10L)
      # The whole curve runs in one background process when run_app() asks for it, since each
      # point the grid adds depends on the power at the last. It is the costlier of the two
      # buttons, and it used to run in this process even under run_app(async = TRUE).
      if (.async_ok) {
        power_result(list(msg = sprintf(
          "Running %d simulations at each of up to 10 sample sizes in a background worker...", n)))
        p <- promises::future_promise({
          power_curve_grid(base_n, cap, n, seeds, function(size, seed) {
            s <- spec; s$units$subject$n <- as.integer(size); s$seed <- seed
            pilotr::power_design(s, n_sims = n)$power
          })
        }, seed = TRUE)
        promises::then(p, onFulfilled = function(curve) report_curve(curve, spec, n),
                          onRejected = power_failed)
      } else {
        # One step of the bar for each point, of the 10 a grid can hold.
        withProgress(message = "Computing the power curve...", value = 0, {
          curve <- power_curve_grid(base_n, cap, n, seeds, function(size, seed) {
            s <- spec; s$units$subject$n <- as.integer(size); s$seed <- seed
            pw <- power_design(s, n_sims = n)$power
            incProgress(1 / length(seeds), detail = sprintf("N = %d", as.integer(size)))
            pw
          })
        })
        report_curve(curve, spec, n)
      }
    }, error = power_failed)
  })

  # Render as text (not print) to keep the message free of a trailing NULL. Shiny skips
  # sending an output whose value is unchanged, so a repeated result (or the async "Running…"
  # placeholder followed by an identical result) could leave a stale line on screen. Toggling
  # an invisible trailing space by the click count makes each value distinct, so the update is
  # always sent; the space does not show in the monospaced output.
  output$power_out <- renderText({
    n <- sum(input$run_power, input$run_curve)
    r <- power_result()
    txt <- if (is.null(r))
      sprintf(paste("Set the number of simulations, then select Run power analysis. The app caps",
                    "the number at %d. For more, call power_design() directly, or set",
                    "PILOTR_MAX_SIMS before launching the app."), MAX_SIMS)
    else if (!is.null(r$msg)) r$msg
    else paste0(
      sprintf("Simulations  : %d\n", r$n_sims),
      sprintf("Power        : %.3f (MCSE %.3f; 95%% Wilson [%.3f, %.3f])\n",
              r$power, r$power_mcse, r$power_lo, r$power_hi),
      "Type S error : ", type_s_text(r), "\n",
      "Type M (exag): ", type_m_text(r), "\n",
      sprintf("True effect  : %.3f | mean estimate: %.3f", r$true_effect, r$mean_estimate))
    paste0(txt, strrep(" ", n %% 2))
  })
  output$power_plot <- renderPlot({
    pc <- power_curve_data(); if (is.null(pc)) return(NULL)
    df <- data.frame(n = pc$grid, power = pc$pw)
    p <- ggplot(df, aes(n, power)) +
      geom_hline(yintercept = 0.8, linetype = 2, colour = "#888888") +
      annotate("text", x = min(df$n), y = 0.8, label = "0.80 target",
               hjust = 0, vjust = -0.6, colour = "#888888", size = 3.6)
    if (!is.null(pc$solved))
      p <- p +
        annotate("rect", xmin = pc$solved$lo, xmax = pc$solved$hi,
                 ymin = -Inf, ymax = Inf, fill = "#888888", alpha = .15) +
        geom_vline(xintercept = pc$solved$value, linetype = 2, colour = "#888888")
    p +
      geom_line(colour = "#2C6FB0", linewidth = 0.9) +
      geom_point(colour = "#2C6FB0", size = 3) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(x = expression(italic(N) ~ "subjects"), y = "Power", title = "Power curve") +
      theme_minimal(base_size = 14)
  })

  output$dl_spec <- downloadHandler(
    filename = function() paste0(current_stem(), ".json"),
    content = function(file) write_text_download(spec_json(current_spec()), file))
  # The data the Data tab shows, named after the specification and seed that made them.
  output$dl_data <- downloadHandler(
    filename = function() {
      s <- sim()$spec
      paste0(file_stem(s$name), "_seed", format(s$seed, scientific = FALSE), ".csv")
    },
    content = function(file) write.csv(sim()$data, file, row.names = FALSE))
}

shinyApp(ui, server)
