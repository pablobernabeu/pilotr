# simdgp no-code app -- the third interface over the shared design spec.
#
# Thin client: every control writes into the portable JSON spec, which downloads and runs
# unchanged in the R and Python packages. Launch with simdgp::run_app() (installed) or
# shiny::runApp("toolkit/r/simdgp/inst/app") (from source).

library(shiny)

# When the package is loaded (run_app), its functions are available; from source, locate
# and source the engine + spec-builder. (Installed packages have no R/ source files, so we
# only source when the functions are not already present.)
if (!exists("simulate_design", mode = "function")) {
  .src <- function(rel, required = TRUE) {
    for (b in c("../../R", "../R", "R", ".")) {
      p <- file.path(b, rel); if (file.exists(p)) { source(p); return(invisible(TRUE)) }
    }
    if (required) stop("cannot find ", rel); invisible(FALSE)
  }
  for (f in c("core.R", "simulate.R", "power.R", "spec_builder.R")) .src(f)
  .src("power_mixed.R", required = FALSE)  # mixed-effects power if available
}

MAX_SIMS <- as.integer(Sys.getenv("SIMDGP_MAX_SIMS", "5000"))
# Async only when running as the installed package with future+promises (workers reload the
# package). From source / serverless this is FALSE and power runs synchronously.
.async_ok <- isNamespaceLoaded("simdgp") &&
  nzchar(system.file(package = "future")) && nzchar(system.file(package = "promises"))

# ---------------------------------------------------------------- UI ----
ui <- fluidPage(
  titlePanel("simdgp — design · simulate · power (one spec, three interfaces)"),
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
                    "Ordinal (Likert, cumulative-logit)" = "ordinal")),
      textInput("resp_name", "Response name (blank = auto)", ""),
      conditionalPanel("input.family == 'gaussian' || input.family == 'shifted_lognormal'",
        numericInput("sigma", "Residual SD (log scale for RT)", 10, min = 0)),
      conditionalPanel("input.family == 'shifted_lognormal'",
        numericInput("shift", "Shift / non-decision time", 200)),
      conditionalPanel("input.family == 'ordinal'",
        textInput("thresholds", "Thresholds (comma-separated)", "-2, -0.6, 0.6, 2")),
      tags$hr(),
      actionButton("simulate", "Simulate", class = "btn-primary"),
      downloadButton("dl_spec", "Download spec (.json)"),
      downloadButton("dl_data", "Download data (.csv)")
    ),
    mainPanel(
      width = 8,
      tabsetPanel(
        tabPanel("Design spec (JSON)",
          p("This portable spec is the single source of truth. Download it and run it ",
            "unchanged in R or Python — you get the identical data set."),
          verbatimTextOutput("json")),
        tabPanel("Data", verbatimTextOutput("dims"), tableOutput("head")),
        tabPanel("Summary & plot", verbatimTextOutput("summary"), plotOutput("plot", height = "320px")),
        tabPanel("Power & design analysis",
          p("Simulation-based power with Type S / Type M (Gelman & Carlin, 2014)."),
          numericInput("n_sims", "Simulations (capped in-app)", 1000, min = 100, max = MAX_SIMS, step = 100),
          actionButton("run_power", "Run power analysis"),
          verbatimTextOutput("power_out")),
        tabPanel("Reproduce in R / Python", verbatimTextOutput("repro"))
      )
    )
  )
)

# ------------------------------------------------------------ server ----
server <- function(input, output, session) {

  current_spec <- reactive({
    build_spec(list(
      name = input$name, seed = input$seed, n_subject = input$n_subject,
      include_items = input$include_items, n_item = input$n_item,
      design_kind = input$design_kind, factor_name = input$factor_name,
      lev1 = input$lev1, lev2 = input$lev2, intercept = input$intercept, effect = input$effect,
      subj_int_sd = input$subj_int_sd, subj_slope_sd = input$subj_slope_sd, subj_corr = input$subj_corr,
      item_int_sd = input$item_int_sd, item_slope_sd = input$item_slope_sd, item_corr = input$item_corr,
      family = input$family, resp_name = input$resp_name, sigma = input$sigma,
      shift = input$shift, thresholds = input$thresholds
    ))
  })

  data <- eventReactive(input$simulate, simulate_design(current_spec()), ignoreNULL = FALSE)

  output$json <- renderText(spec_json(current_spec()))
  output$dims <- renderText({ d <- data(); sprintf("Simulated %d rows x %d columns (seed %d).", nrow(d), ncol(d), input$seed) })
  output$head <- renderTable(head(data(), 10))

  output$summary <- renderPrint({
    d <- data(); yn <- current_spec()$response$name; fn <- input$factor_name
    if (is.numeric(d[[yn]])) {
      agg <- aggregate(d[[yn]], list(d[[fn]]), function(x) c(mean = mean(x), sd = sd(x), n = length(x)))
      cat("Mean (SD) of", yn, "by", fn, ":\n"); print(do.call(data.frame, agg))
    } else { cat("Counts of", yn, "by", fn, ":\n"); print(table(d[[fn]], d[[yn]])) }
  })

  output$plot <- renderPlot({
    d <- data(); yn <- current_spec()$response$name; fn <- input$factor_name
    if (is.numeric(d[[yn]])) boxplot(d[[yn]] ~ d[[fn]], xlab = fn, ylab = yn,
                                     col = c("#2C6FB0", "#B0402C"), main = paste("Distribution of", yn))
    else barplot(table(d[[fn]], d[[yn]]), beside = TRUE, legend = TRUE,
                 col = c("#2C6FB0", "#B0402C"), main = paste("Counts of", yn))
  })

  output$repro <- renderText(paste0(
    "# Download the spec (first tab) as design.json, then:\n\n",
    "# --- R ---\nlibrary(simdgp)\nd <- simulate_design(\"design.json\")\n\n",
    "# --- Python ---\nfrom simdgp import simulate\nd = simulate(\"design.json\")\nd.to_csv(\"data.csv\")\n\n",
    "# Same spec + seed => identical data in all three interfaces."))

  # ---- power: capped, async when installed (worker process), else synchronous ----
  power_result <- reactiveVal(NULL)
  observeEvent(input$run_power, {
    spec <- current_spec()
    if (spec$response$family != "gaussian" || identical(input$design_kind, "within")) {
      power_result(list(msg = paste0("The in-app power backend covers the two-group Gaussian design.\n",
        "For crossed mixed-effects designs use power_mixed() in the package (fits lme4; minutes).")))
      return()
    }
    n <- min(max(as.integer(input$n_sims), 100L), MAX_SIMS)
    if (.async_ok) {
      power_result(list(msg = sprintf("Running %d simulations in a background worker...", n)))
      p <- promises::future_promise({ simdgp::power_design(spec, n_sims = n) }, seed = TRUE)
      promises::then(p, onFulfilled = function(res) power_result(res),
                        onRejected = function(e) power_result(list(msg = paste("error:", conditionMessage(e)))))
    } else {
      withProgress(message = sprintf("Simulating %d datasets...", n), value = 0.5,
                   power_result(power_design(spec, n_sims = n)))
    }
  })
  output$power_out <- renderPrint({
    r <- power_result()
    if (is.null(r)) { cat(sprintf("Set simulations (capped at %d here; unlimited when you install the package) and click Run.", MAX_SIMS)); return() }
    if (!is.null(r$msg)) { cat(r$msg); return() }
    cat(sprintf("Simulations  : %d\nPower        : %.3f\nType S error : %.4f\nType M (exag): %.3f\nTrue effect  : %.3f | mean estimate: %.3f\n",
                r$n_sims, r$power, r$type_s, r$type_m, r$true_effect, r$mean_estimate))
  })

  output$dl_spec <- downloadHandler(
    filename = function() paste0(input$name, ".json"),
    content = function(file) writeLines(spec_json(current_spec()), file))
  output$dl_data <- downloadHandler(
    filename = function() paste0(input$name, "_seed", input$seed, ".csv"),
    content = function(file) write.csv(simulate_design(current_spec()), file, row.names = FALSE))
}

shinyApp(ui, server)
