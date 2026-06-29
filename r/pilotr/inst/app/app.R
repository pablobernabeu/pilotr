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
  .resolved <- character(0)
  .src <- function(rel, required = TRUE) {
    for (b in c("../../R", "../R", "R", ".")) {
      p <- file.path(b, rel); if (file.exists(p)) { source(p); .resolved <<- c(.resolved, normalizePath(p)); return(invisible(TRUE)) }
    }
    if (required) stop("cannot find ", rel); invisible(FALSE)
  }
  for (f in c("core.R", "simulate.R", "power.R", "spec_builder.R")) .src(f)
  .src("power_mixed.R", required = FALSE)  # mixed-effects power if available
  ENGINE_FILES <- .resolved
}

MAX_SIMS <- as.integer(Sys.getenv("PILOTR_MAX_SIMS", "5000"))
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
           "sample size. Crossed mixed-effects power (via lme4) and precision/ROPE design ",
           "analysis run in the installed package; the specification you build here drives all ",
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
        textAreaInput("spec_json_in", NULL, "", rows = 5,
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
          numericInput("n_sims", "Simulations (capped in-app)", 1000, min = 100, max = MAX_SIMS, step = 100),
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

  current_spec <- reactive({
    txt <- input$spec_json_in
    if (!is.null(txt) && nzchar(trimws(txt))) {                 # advanced override
      parsed <- tryCatch(jsonlite::fromJSON(txt, simplifyVector = TRUE,
                          simplifyDataFrame = FALSE, simplifyMatrix = FALSE), error = function(e) NULL)
      if (!is.null(parsed) && !is.null(parsed$response)) return(parsed)
    }
    build_spec(list(
      name = input$name, seed = input$seed, n_subject = input$n_subject,
      include_items = input$include_items, n_item = input$n_item,
      design_kind = input$design_kind, factor_name = input$factor_name,
      lev1 = input$lev1, lev2 = input$lev2, intercept = input$intercept, effect = input$effect,
      subj_int_sd = input$subj_int_sd, subj_slope_sd = input$subj_slope_sd, subj_corr = input$subj_corr,
      item_int_sd = input$item_int_sd, item_slope_sd = input$item_slope_sd, item_corr = input$item_corr,
      family = input$family, resp_name = input$resp_name, sigma = input$sigma,
      shift = input$shift, thresholds = input$thresholds, phi = input$phi
    ))
  })

  data <- eventReactive(input$simulate, simulate_design(current_spec()), ignoreNULL = FALSE)

  output$json <- renderText(spec_json(current_spec()))
  output$dims <- renderText({ d <- data(); sprintf("Simulated %d rows x %d columns (seed %d).", nrow(d), ncol(d), input$seed) })
  output$head <- renderTable(head(data(), 10), striped = TRUE, spacing = "xs")

  output$summary <- renderPrint({
    d <- data(); spec <- current_spec(); yn <- spec$response$name; fn <- .grp_col(spec)
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
    d <- data(); spec <- current_spec(); yn <- spec$response$name; fn <- .grp_col(spec)
    has_grp <- !is.null(fn) && fn %in% names(d)
    pal2 <- c("#2C6FB0", "#B0402C")
    base <- theme_minimal(base_size = 14)
    if (is.numeric(d[[yn]])) {
      if (has_grp)
        ggplot(d, aes(.data[[fn]], .data[[yn]], fill = .data[[fn]])) +
          geom_boxplot(alpha = 0.85, outlier.alpha = 0.35) +
          scale_fill_manual(values = pal2, guide = "none") +
          labs(x = fn, y = yn, title = paste("Distribution of", yn)) + base
      else
        ggplot(d, aes(.data[[yn]])) +
          geom_histogram(bins = 30, fill = "#2C6FB0", colour = "white") +
          labs(x = yn, y = "count", title = paste("Distribution of", yn)) + base
    } else if (has_grp)
      ggplot(d, aes(.data[[yn]], fill = .data[[fn]])) +
        geom_bar(position = "dodge") +
        scale_fill_manual(values = pal2, name = fn) +
        labs(x = yn, y = "count", title = paste("Counts of", yn)) + base
    else
      ggplot(d, aes(.data[[yn]])) +
        geom_bar(fill = "#2C6FB0") +
        labs(x = yn, y = "count", title = paste("Counts of", yn)) + base
  })

  output$rscript <- renderText(generate_r_script(current_spec()))
  output$repro_py <- renderText(paste0(
    "# The same design also runs in Python (bit-identical given the same seed):\n",
    "from pilotr import simulate\n",
    "d = simulate(\"design.json\")   # download the spec from the first tab\n",
    "d.to_csv(\"data.csv\")"))

  output$dl_rscript <- downloadHandler(
    filename = function() paste0(input$name, ".R"),
    content = function(file) writeLines(generate_r_script(current_spec()), file))

  # ---- Verify: run the exported design in a clean R subprocess and compare ----
  verify_result <- reactiveVal(NULL)
  observeEvent(input$verify_code, {
    spec <- current_spec()
    ref <- simulate_design(spec); yn <- spec$response$name
    ref_chk <- if (is.numeric(ref[[yn]])) sum(ref[[yn]]) else paste(ref[[yn]], collapse = "")
    if (!requireNamespace("callr", quietly = TRUE)) {
      verify_result(list(msg = "Install the 'callr' package to verify in a clean R session.")); return()
    }
    withProgress(message = "Running the design in a clean R session...", value = 0.5, {
      res <- tryCatch(callr::r(function(json, files) {
        if (is.null(files)) library(pilotr) else for (f in files) source(f)
        s <- jsonlite::fromJSON(json, simplifyVector = TRUE, simplifyDataFrame = FALSE, simplifyMatrix = FALSE)
        d <- simulate_design(s); yn <- s$response$name
        list(n = nrow(d), chk = if (is.numeric(d[[yn]])) sum(d[[yn]]) else paste(d[[yn]], collapse = ""))
      }, args = list(json = spec_json(spec), files = ENGINE_FILES)), error = function(e) e)
    })
    if (inherits(res, "error")) { verify_result(list(msg = paste("error:", conditionMessage(res)))); return() }
    ident <- isTRUE(all.equal(res$chk, ref_chk)) && res$n == nrow(ref)
    verify_result(list(ok = ident, n = res$n, ref_n = nrow(ref)))
  })
  output$verify_out <- renderText({
    r <- verify_result()
    if (is.null(r)) return("Select Verify to run the script in a fresh R process and confirm that it reproduces this data.")
    if (!is.null(r$msg)) return(r$msg)
    if (isTRUE(r$ok)) sprintf("Reproduces identically in a clean R session.\n  %d rows (app) match %d rows (clean run), and the response checksum matches.", r$ref_n, r$n)
    else sprintf("Mismatch. The app produced %d rows and the clean run produced %d rows, or the checksum differs.", r$ref_n, r$n)
  })

  # ---- power: point estimate + curve, capped, async when installed (worker process) ----
  power_result     <- reactiveVal(NULL)
  power_curve_data <- reactiveVal(NULL)
  gaussian_two_group <- function(spec)
    identical(spec$response$family, "gaussian") &&
      length(spec$factors) >= 1 && !is.null(spec$factors[[1]]$between)
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
    if (!gaussian_two_group(spec)) { power_result(list(msg = not_supported)); return() }
    n <- min(max(as.integer(input$n_sims), 100L), MAX_SIMS)
    if (.async_ok) {
      power_result(list(msg = sprintf("Running %d simulations in a background worker...", n)))
      p <- promises::future_promise({ pilotr::power_design(spec, n_sims = n) }, seed = TRUE)
      promises::then(p, onFulfilled = function(res) power_result(res),
                        onRejected = function(e) power_result(list(msg = paste("error:", conditionMessage(e)))))
    } else {
      withProgress(message = sprintf("Simulating %d datasets...", n), value = 0.5,
                   power_result(power_design(spec, n_sims = n)))
    }
  })

  observeEvent(input$run_curve, {
    spec <- current_spec()
    if (!gaussian_two_group(spec)) { power_result(list(msg = not_supported)); power_curve_data(NULL); return() }
    n <- min(max(as.integer(input$n_sims), 100L), MAX_SIMS)
    base_n <- spec$units$subject$n
    grid <- unique(round(base_n * c(0.5, 0.75, 1, 1.5, 2))); grid <- grid[grid >= 4]
    withProgress(message = "Computing the power curve...", value = 0.3, {
      pw <- vapply(grid, function(nn) {
        s <- spec; s$units$subject$n <- as.integer(nn); power_design(s, n_sims = n)$power
      }, numeric(1))
    })
    power_result(list(msg = sprintf("Power curve at n_sims = %d per point. N subjects = %s.",
                                    n, paste(grid, collapse = ", "))))
    power_curve_data(list(grid = grid, pw = pw))
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
      sprintf("Set the number of simulations (capped at %d in the app, unlimited once you install the package), then select Run power analysis.", MAX_SIMS)
    else if (!is.null(r$msg)) r$msg
    else sprintf("Simulations  : %d\nPower        : %.3f\nType S error : %.4f\nType M (exag): %.3f\nTrue effect  : %.3f | mean estimate: %.3f",
                 r$n_sims, r$power, r$type_s, r$type_m, r$true_effect, r$mean_estimate)
    paste0(txt, strrep(" ", n %% 2))
  })
  output$power_plot <- renderPlot({
    pc <- power_curve_data(); if (is.null(pc)) return(NULL)
    df <- data.frame(n = pc$grid, power = pc$pw)
    ggplot(df, aes(n, power)) +
      geom_hline(yintercept = 0.8, linetype = 2, colour = "#888888") +
      annotate("text", x = min(df$n), y = 0.8, label = "0.80 target",
               hjust = 0, vjust = -0.6, colour = "#888888", size = 3.6) +
      geom_line(colour = "#2C6FB0", linewidth = 0.9) +
      geom_point(colour = "#2C6FB0", size = 3) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(x = expression(italic(N) ~ "subjects"), y = "Power", title = "Power curve") +
      theme_minimal(base_size = 14)
  })

  output$dl_spec <- downloadHandler(
    filename = function() paste0(input$name, ".json"),
    content = function(file) writeLines(spec_json(current_spec()), file))
  output$dl_data <- downloadHandler(
    filename = function() paste0(input$name, "_seed", input$seed, ".csv"),
    content = function(file) write.csv(simulate_design(current_spec()), file, row.names = FALSE))
}

shinyApp(ui, server)
