# pilotr no-code app. This is the third interface over the shared design spec.
#
# A thin client. Every control writes into the portable JSON spec, which downloads and runs
# unchanged in the R and Python packages. Launch with pilotr::run_app() (installed), or with
# shiny::runApp("r/pilotr/inst/app") (from source).

library(shiny)
library(ggplot2)

# When the package is loaded (run_app), its functions are available; from source, locate
# and source the engine + spec-builder. (Installed packages have no R/ source files, so we
# only source when the functions are not already present.)
ENGINE_FILES <- NULL   # resolved source paths (dev) so the verifier can rebuild in a clean R session
if (!exists("simulate_design", mode = "function")) {
  # Load the whole engine directory rather than naming its files. The list this
  # replaced had to track R/ by hand, and when validate.R arrived with 0.3.0 it
  # did not: .SPEC_VERSION went undefined and every Simulate click failed. Each
  # file in R/ only defines functions and constants, so order does not matter.
  # "." last, for a layout that stages the engine beside the app; app.R and the
  # roxygen sentinel are dropped so that case cannot source this file into itself.
  .dir <- Find(dir.exists, c("../../R", "../R", "R", "."))
  if (is.null(.dir)) stop("cannot find the pilotr engine sources")
  .files <- sort(list.files(.dir, pattern = "[.][Rr]$", full.names = TRUE))
  .files <- .files[!basename(.files) %in% c("pilotr-package.R", "app.R")]
  if (!length(.files)) stop("no pilotr engine sources in ", .dir)
  for (f in .files) source(f)
  ENGINE_FILES <- normalizePath(.files)
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (is.character(a) && !nzchar(a))) b else a

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

# Async only when running as the installed package with future+promises (workers reload the
# package). From source / serverless this is FALSE and power runs synchronously.
.async_ok <- isNamespaceLoaded("pilotr") &&
  nzchar(system.file(package = "future")) && nzchar(system.file(package = "promises"))

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

# Fill colours for the plots. A pasted spec may carry a factor with any number of levels, so
# the scales size the palette to the data: a fixed pair handed to a three-level scale left the
# Summary & plot tab showing a ggplot2 error in place of the plot.
PALETTE <- c("#2C6FB0", "#B0402C", "#2E8B57", "#8E6FB0", "#C8922A", "#5A5A5A")
pal <- function(k) if (k <= length(PALETTE)) PALETTE[seq_len(k)] else grDevices::hcl.colors(k, "Dynamic")

# The same theme the lite app sets, so the two variants and the documentation site share one
# identity. A bare fluidPage() renders under Shiny's default Bootstrap 3, whose primary is a
# different blue and which does not define the `mb-3` spacing set below the power buttons.
# Keep this call and app-lite's in step. Without bslib the app still runs, unthemed.
APP_THEME <- if (requireNamespace("bslib", quietly = TRUE))
  bslib::bs_theme(version = 5, primary = "#2C6FB0", "border-radius" = "0.5rem") else NULL

guide_tab <- tabPanel(
  "Guide",
  tags$div(
    style = "max-width: 52rem;",
    tags$p(tags$b("pilotr"), " builds a portable design specification and simulates ",
           "experimental and behavioural data from it. This is the no-code interface. The ",
           "same specification runs unchanged in the R and Python packages and reproduces ",
           "the identical data."),
    tags$ol(
      tags$li("Describe the design on the left: the sample sizes, the factor and its two ",
              "levels, the fixed intercept and effect, and a response family. Changing the ",
              "family resets the intercept and effect to sensible values for that family's ",
              "scale. For within-subjects and crossed designs, add by-subject and by-item ",
              "random effects."),
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
    tags$p("The in-app backend estimates power for the two-group Gaussian design, reports the ",
           "Type S and Type M errors of Gelman and Carlin (2014), and draws a power curve over ",
           "sample size. It then solves that curve for the sample size at which power reaches ",
           "0.80 and reports a confidence interval on it. The heavier analyses stay with the ",
           "packages themselves: crossed mixed-effects power (via lme4) in R and in Python, and ",
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
  theme = APP_THEME,
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
        conditionalPanel("input.include_items", numericInput("n_item", "N items", 24, min = 2))
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
        # A label, not NULL: an empty one names the control for nobody, and the summary above
        # is a sibling rather than a label, so a screen reader had only the placeholder to read.
        textAreaInput("spec_json_in", "Design spec (JSON)", "", rows = 5,
          placeholder = "Paste a pilotr spec with continuous predictors / interactions (e.g. a reading-time design)")),
      actionButton("simulate", "Simulate", class = "btn-primary"),
      downloadButton("dl_spec", "Download spec (.json)"),
      downloadButton("dl_data", "Download data (.csv)")
    ),
    mainPanel(
      width = 8,
      tabsetPanel(
        id = "tabs",
        guide_tab,
        tabPanel("Design spec (JSON)",
          p("This portable spec is the single source of truth. Download it and run it ",
            "unchanged in R or Python to obtain the identical data set."),
          verbatimTextOutput("json")),
        tabPanel("Data", verbatimTextOutput("dims"), tableOutput("head")),
        tabPanel("Summary & plot", verbatimTextOutput("summary"), plotOutput("plot", height = "320px")),
        tabPanel("Power & design analysis",
          p("Simulation-based power with Type S / Type M (Gelman & Carlin, 2014), for the ",
            "two-group Gaussian design."),
          numericInput("n_sims", "Simulations (capped in-app)", N_SIMS_DEFAULT,
                       min = N_SIMS_MIN, max = MAX_SIMS, step = 100),
          div(class = "mb-3",
              actionButton("run_power", "Run power analysis", class = "btn-primary"), " ",
              actionButton("run_curve", "Power curve over N")),
          verbatimTextOutput("power_out"),
          plotOutput("power_plot", height = "300px")),
        tabPanel("Reproducible R script",
          p("Your no-code design as a self-contained R script. You can download it, or ",
            "verify that it reproduces this exact data by running it in a clean R session."),
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

  parse_error <- reactiveVal(NULL)

  # Returns the spec, or NULL if a pasted spec is invalid, with the message in parse_error().
  # The message travels in a reactive value rather than a validate() call so that the observers
  # and the downloads can report it too, as the lite app's already did: validate() inside an
  # observer aborts it silently, and inside a downloadHandler it fails the browser's request.
  current_spec <- reactive({
    txt <- input$spec_json_in
    if (!is.null(txt) && nzchar(trimws(txt))) {                 # advanced override
      # As in the lite app, route the paste through load_spec() (staged in a temp file) rather
      # than a bare fromJSON, so a pasted spec faces the same validation as a loaded one.
      # Several of the ways a spec can be wrong produce plausible data rather than an error
      # (a mistyped coefficient key resolves to no column and silently generates null-effect
      # data), which is exactly what 0.3.0's validation refuses.
      spec <- tryCatch({
        tf <- tempfile(fileext = ".json"); writeLines(txt, tf); load_spec(tf)
      }, error = function(e) e)
      if (inherits(spec, "error")) { parse_error(conditionMessage(spec)); return(NULL) }
      parse_error(NULL)
      return(spec)
    }
    # The controls can also build a specification the package refuses: a cleared numeric box
    # arrives as NA, and a response named after the factor takes its column. Validating the
    # built spec here is what keeps the Design spec tab, which the app calls the single source
    # of truth, from showing and downloading a file that load_spec() will not read.
    spec <- tryCatch(validate_spec(build_spec(list(
      name = input$name, seed = input$seed, n_subject = input$n_subject,
      include_items = input$include_items, n_item = input$n_item,
      design_kind = input$design_kind, factor_name = input$factor_name,
      lev1 = input$lev1, lev2 = input$lev2, intercept = input$intercept, effect = input$effect,
      subj_int_sd = input$subj_int_sd, subj_slope_sd = input$subj_slope_sd, subj_corr = input$subj_corr,
      item_int_sd = input$item_int_sd, item_slope_sd = input$item_slope_sd, item_corr = input$item_corr,
      family = input$family, resp_name = input$resp_name, sigma = input$sigma,
      shift = input$shift, thresholds = input$thresholds, phi = input$phi
    ))), error = function(e) e)
    if (inherits(spec, "error")) { parse_error(conditionMessage(spec)); return(NULL) }
    parse_error(NULL)
    spec
  })

  spec_req <- function() {
    s <- current_spec()
    validate(need(!is.null(s), parse_error() %||% "Please enter a valid design specification."))
    s
  }

  # Simulate captures the specification it ran alongside the data, so that every tab describes
  # one snapshot. Reading the response and factor names from the live specification instead put
  # the names out of step with the columns as soon as a control was changed without a re-click,
  # and the summary and the plot then failed with a raw R message.
  #
  # A specification the engine refuses, such as one whose response repeats the factor name, comes
  # back as a message rather than an error, so that it reaches the tabs as a sentence.
  snapshot <- eventReactive(input$simulate, {
    s <- current_spec(); if (is.null(s)) return(NULL)
    tryCatch(list(spec = s, data = simulate_design(s)),
             error = function(e) list(error = conditionMessage(e)))
  }, ignoreNULL = FALSE)

  snap_req <- function() {
    snap <- snapshot()
    validate(need(!is.null(snap) && is.null(snap$error),
                  snap$error %||% parse_error() %||%
                    "Please correct the specification, then select Simulate."))
    snap
  }
  data_req <- function() snap_req()$data

  output$json <- renderText(spec_json(spec_req()))
  # The notice is what keeps the Data tab honest once the design has moved on: the table below
  # still shows the snapshot, not the specification the other tabs describe.
  output$dims <- renderText({
    snap <- snap_req(); d <- snap$data; live <- current_spec()
    sprintf("Simulated %d rows x %d columns (seed %s).%s", nrow(d), ncol(d), snap$spec$seed,
            if (!is.null(live) && !identical(live, snap$spec))
              " The design has changed since then; select Simulate to bring this up to date." else "")
  })
  output$head <- renderTable(head(data_req(), 10), striped = TRUE, spacing = "xs")

  output$summary <- renderPrint({
    snap <- snap_req(); d <- snap$data; spec <- snap$spec
    yn <- spec$response$name; fn <- .grp_col(spec)
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
    snap <- snap_req(); d <- snap$data; spec <- snap$spec
    yn <- spec$response$name; fn <- .grp_col(spec)
    has_grp <- !is.null(fn) && fn %in% names(d)
    base <- theme_minimal(base_size = 14)
    if (is.numeric(d[[yn]])) {
      if (has_grp)
        ggplot(d, aes(.data[[fn]], .data[[yn]], fill = .data[[fn]])) +
          geom_boxplot(alpha = 0.85, outlier.alpha = 0.35) +
          scale_fill_manual(values = pal(nlevels(factor(d[[fn]]))), guide = "none") +
          labs(x = fn, y = yn, title = paste("Distribution of", yn)) + base
      else
        ggplot(d, aes(.data[[yn]])) +
          geom_histogram(bins = 30, fill = PALETTE[1], colour = "white") +
          labs(x = yn, y = "count", title = paste("Distribution of", yn)) + base
    } else if (has_grp)
      ggplot(d, aes(.data[[yn]], fill = .data[[fn]])) +
        geom_bar(position = "dodge") +
        scale_fill_manual(values = pal(nlevels(factor(d[[fn]]))), name = fn) +
        labs(x = yn, y = "count", title = paste("Counts of", yn)) + base
    else
      ggplot(d, aes(.data[[yn]])) +
        geom_bar(fill = PALETTE[1]) +
        labs(x = yn, y = "count", title = paste("Counts of", yn)) + base
  })

  output$rscript <- renderText(generate_r_script(spec_req()))
  output$repro_py <- renderText(paste0(
    "# The same design also runs in Python (bit-identical given the same seed):\n",
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

  # Downloads are named after the specification they carry rather than after the sidebar boxes,
  # as in the lite app. With a pasted spec in force the two are different designs, and a data
  # file named after a seed that did not generate it is worse than one carrying no seed at all,
  # the filename being all that travels with a CSV once it leaves the app. An emptied Design
  # name falls back to 'design' rather than writing a dot-file.
  download_name <- function(ext, seed = FALSE) {
    s <- current_spec()
    paste0(s$name %||% "design",
           if (seed && !is.null(s$seed)) paste0("_seed", s$seed) else "", ext)
  }

  # The three downloads degrade gracefully on an invalid pasted spec, as in the lite app: a
  # downloadHandler that stops leaves the browser with a failed request and no explanation.
  output$dl_rscript <- downloadHandler(
    filename = function() download_name(".R"),
    content = function(file) {
      s <- current_spec()
      oops <- "# Invalid specification. Correct it in the app, then download the script again."
      write_text_download(if (is.null(s)) oops else generate_r_script(s), file)
    })

  # ---- Verify: run the exported R script in a clean R subprocess and compare ----
  verify_result <- reactiveVal(NULL)
  observeEvent(input$verify_code, {
    spec <- current_spec()
    if (is.null(spec)) {
      verify_result(list(msg = parse_error() %||% "Please correct the specification first.")); return()
    }
    ref <- tryCatch(simulate_design(spec), error = function(e) e)
    if (inherits(ref, "error")) { verify_result(list(msg = conditionMessage(ref))); return() }
    if (!requireNamespace("callr", quietly = TRUE)) {
      verify_result(list(msg = "Install the 'callr' package to verify in a clean R session.")); return()
    }
    # The subprocess runs the script this tab offers, which is what the button promises. Handing
    # it the JSON instead exercised the JSON writer and left the script's own writer, which
    # quotes list names and formats every number itself, unverified: a defect confined to it
    # would have left the button green. The comparison is on the whole data frame, so
    # 'reproduces this data' covers the column names, the types and the row order too.
    withProgress(message = "Running the exported script in a clean R session...", value = 0.5, {
      res <- tryCatch(callr::r(function(script, files) {
        env <- new.env(parent = globalenv())
        # From source there is no installed package to attach, and the engine files the app
        # itself read are the code the script has to reproduce with, so its library() call
        # stands down once they are in place.
        if (is.null(files)) library(pilotr)
        else { for (f in files) source(f, local = env); env$library <- function(...) invisible(NULL) }
        eval(parse(text = script), envir = env)
        env$data
      }, args = list(script = generate_r_script(spec), files = ENGINE_FILES)), error = function(e) e)
    })
    if (inherits(res, "error")) { verify_result(list(msg = paste("error:", conditionMessage(res)))); return() }
    if (!is.data.frame(res)) { verify_result(list(msg = "The script ran but left no data set behind.")); return() }
    verify_result(list(ok = identical(res, ref), n = nrow(res), ref_n = nrow(ref)))
  })
  output$verify_out <- renderText({
    r <- verify_result()
    if (is.null(r)) return("Select Verify to run the script in a fresh R process and confirm that it reproduces this data.")
    if (!is.null(r$msg)) return(r$msg)
    if (isTRUE(r$ok)) sprintf("The script reproduces this data bit for bit in a clean R session.\n  %d rows, every column identical.", r$ref_n)
    else sprintf("Mismatch. The app produced %d rows and the clean run of the script produced %d; the two data sets are not identical.", r$ref_n, r$n)
  })

  # ---- power: point estimate + curve, capped, async when installed (worker process) ----
  power_result     <- reactiveVal(NULL)
  power_curve_data <- reactiveVal(NULL)
  # The shape power_design() covers: a gaussian response and exactly one 2-level between
  # factor. Testing only the first factor let a 3-level design (or a second between factor)
  # through to the engine, whose refusal then travelled as an unhandled error.
  gaussian_two_group <- function(spec) {
    between <- Filter(function(f) !is.null(f$between), spec$factors)
    identical(spec$response$family, "gaussian") &&
      length(between) == 1L && length(between[[1]]$levels) == 2L
  }
  not_supported <- paste0(
    "The in-app power backend covers the two-group Gaussian design. For a crossed\n",
    "mixed-effects design, download the spec (the Design spec tab) and run it directly:\n\n",
    "R (lme4; may take a few minutes):\n",
    "    library(pilotr)\n",
    "    spec <- load_spec(\"design.json\")\n",
    "    power_mixed(spec, n_sims = 200)\n",
    "    power_curve_mixed(spec, subject_ns = c(20, 40, 60), n_sims = 200)\n\n",
    "Python (statsmodels backend):\n",
    "    from pilotr import load_spec, power_mixed\n",
    "    power_mixed(load_spec(\"design.json\"), n_sims=200)")

  observeEvent(input$run_power, {
    power_curve_data(NULL)
    spec <- current_spec()
    if (is.null(spec)) { power_result(list(msg = parse_error() %||% "Please correct the specification first.")); return() }
    if (!gaussian_two_group(spec)) { power_result(list(msg = not_supported)); return() }
    n <- .n_sims_input(input$n_sims)
    if (.async_ok) {
      power_result(list(msg = sprintf("Running %d simulations in a background worker...", n)))
      p <- promises::future_promise({ pilotr::power_design(spec, n_sims = n) }, seed = TRUE)
      promises::then(p, onFulfilled = function(res) power_result(res),
                        onRejected = function(e) power_result(list(msg = paste("error:", conditionMessage(e)))))
    } else {
      # An error in an observer ends the session, and with it the design the user built, so
      # every engine refusal that the guard above does not cover (too few subjects per group,
      # say) is reported as a line of text instead.
      withProgress(message = sprintf("Simulating %d datasets...", n), value = 0.5,
                   power_result(tryCatch(power_design(spec, n_sims = n),
                                         error = function(e) list(msg = paste("error:", conditionMessage(e))))))
    }
  })

  observeEvent(input$run_curve, {
    spec <- current_spec()
    if (is.null(spec)) { power_result(list(msg = parse_error() %||% "Please correct the specification first.")); power_curve_data(NULL); return() }
    if (!gaussian_two_group(spec)) { power_result(list(msg = not_supported)); power_curve_data(NULL); return() }
    n <- .n_sims_input(input$n_sims)
    base_n <- spec$units$subject$n
    grid <- unique(round(base_n * c(0.5, 0.75, 1, 1.5, 2))); grid <- grid[grid >= 4]
    if (!length(grid)) {
      power_result(list(msg = "The power curve needs at least 4 subjects. Raise N subjects, then select Power curve."))
      power_curve_data(NULL); return()
    }
    # The curve is the sweep the point estimate repeats once per grid point, so it is the
    # more expensive of the two buttons and the one that most needs the background worker
    # run_app(async = TRUE) promises. Reporting the solved curve is shared by both branches.
    report <- function(pw) {
      if (inherits(pw, "error")) {
        power_result(list(msg = paste("error:", conditionMessage(pw)))); power_curve_data(NULL); return()
      }
      # A dashed target line leaves the reader to judge the crossing. target_n() estimates
      # it instead, and reports its refusal when the curve does not settle the question,
      # which is more use than a number the sweep cannot support.
      solved <- tryCatch(target_n(data.frame(n_subject = grid, power = pw, n_sims = n),
                                  target = 0.8),
                         error = function(e) conditionMessage(e))
      power_result(list(msg = paste0(
        sprintf("Power curve at n_sims = %d per point. N subjects = %s.\n",
                n, paste(grid, collapse = ", ")),
        if (is.character(solved)) paste0("Target power 0.80: ", solved)
        else sprintf("Target power 0.80: N = %d subjects (95%% interval %d to %d).",
                     solved$n, solved$n_lo, solved$n_hi))))
      power_curve_data(list(grid = grid, pw = pw,
                            solved = if (is.character(solved)) NULL else solved))
    }
    if (.async_ok) {
      power_curve_data(NULL)
      power_result(list(msg = sprintf("Running %d simulations at each of %d sample sizes in a background worker...",
                                      n, length(grid))))
      p <- promises::future_promise({
        vapply(grid, function(nn) {
          s <- spec; s$units$subject$n <- as.integer(nn); pilotr::power_design(s, n_sims = n)$power
        }, numeric(1))
      }, seed = TRUE)
      promises::then(p, onFulfilled = report,
                        onRejected = function(e) {
                          power_result(list(msg = paste("error:", conditionMessage(e))))
                          power_curve_data(NULL)
                        })
    } else {
      # One step per grid point rather than one for the whole sweep, which left the bar at a
      # third of the way across for as long as the sweep took.
      report(withProgress(message = "Computing the power curve...", value = 0,
        tryCatch(vapply(seq_along(grid), function(i) {
          s <- spec; s$units$subject$n <- as.integer(grid[i])
          pw_i <- power_design(s, n_sims = n)$power
          incProgress(1 / length(grid),
                      detail = sprintf("N = %d (%d of %d)", grid[i], i, length(grid)))
          pw_i
        }, numeric(1)), error = function(e) e)))
    }
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
      sprintf("Set the number of simulations (capped at %d here; power_design() takes any count when called directly), then select Run power analysis.", MAX_SIMS)
    else if (!is.null(r$msg)) r$msg
    else sprintf("Simulations  : %d\nPower        : %.3f\nType S error : %.4f\nType M (exag): %.3f\nTrue effect  : %.3f | mean estimate: %.3f",
                 r$n_sims, r$power, r$type_s, r$type_m, r$true_effect, r$mean_estimate)
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
      geom_line(colour = PALETTE[1], linewidth = 0.9) +
      geom_point(colour = PALETTE[1], size = 3) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(x = expression(italic(N) ~ "subjects"), y = "Power", title = "Power curve") +
      theme_minimal(base_size = 14)
  })

  output$dl_spec <- downloadHandler(
    filename = function() download_name(".json"),
    content = function(file) {
      s <- current_spec()
      write_text_download(if (is.null(s)) "{}" else spec_json(s), file)
    })
  output$dl_data <- downloadHandler(
    filename = function() download_name(".csv", seed = TRUE),
    content = function(file) {
      s <- current_spec()
      d <- if (is.null(s)) NULL else tryCatch(simulate_design(s), error = function(e) NULL)
      if (is.null(d)) writeLines("", file) else write.csv(d, file, row.names = FALSE)
    })
}

shinyApp(ui, server)
