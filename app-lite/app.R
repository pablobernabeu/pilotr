# simdgp LITE app -- the serverless (shinylive / webR) build of the light path.
#
# Runs entirely in the browser via WebAssembly: zero server, infinite concurrent users.
# Covers design building, simulation, summaries, the Gaussian power backend, and spec/data
# download. Heavy crossed mixed-effects power (lme4) is intentionally a local-install
# feature. The engine + spec-builder files are staged alongside this app by build.R.

library(shiny)
for (f in c("core.R", "simulate.R", "power.R", "spec_builder.R")) source(f)

ui <- fluidPage(
  titlePanel("simdgp (lite) — design & simulate in your browser"),
  tags$p(tags$em("Runs fully in your browser (no server). For crossed mixed-effects power, ",
                 "install the simdgp R or Python package.")),
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
        checkboxInput("include_items", "Crossed with items", TRUE),
        conditionalPanel("input.include_items", numericInput("n_item", "N items", 24, min = 2)),
        tags$b("By-subject random effects"),
        fluidRow(
          column(4, numericInput("subj_int_sd", "Int SD", 0.12, min = 0)),
          column(4, numericInput("subj_slope_sd", "Slope SD", 0.04, min = 0)),
          column(4, numericInput("subj_corr", "corr", 0.2, min = -1, max = 1))
        ),
        conditionalPanel("input.include_items",
          tags$b("By-item random effects"),
          fluidRow(
            column(4, numericInput("item_int_sd", "Int SD", 0.08, min = 0)),
            column(4, numericInput("item_slope_sd", "Slope SD", 0.02, min = 0)),
            column(4, numericInput("item_corr", "corr", -0.1, min = -1, max = 1))
          ))
      ),
      tags$hr(),
      fluidRow(
        column(6, textInput("factor_name", "Factor", "group")),
        column(3, textInput("lev1", "Level 1", "control")),
        column(3, textInput("lev2", "Level 2", "treatment"))
      ),
      fluidRow(
        column(6, numericInput("intercept", "Intercept", 100)),
        column(6, numericInput("effect", "Effect (on -0.5/+0.5)", 5))
      ),
      selectInput("family", "Response family",
                  c("Gaussian" = "gaussian", "Shifted lognormal (RT)" = "shifted_lognormal",
                    "Bernoulli (logit)" = "bernoulli", "Poisson (log)" = "poisson",
                    "Ordinal (Likert)" = "ordinal")),
      textInput("resp_name", "Response name (blank = auto)", ""),
      conditionalPanel("input.family == 'gaussian' || input.family == 'shifted_lognormal'",
        numericInput("sigma", "Residual SD", 10, min = 0)),
      conditionalPanel("input.family == 'shifted_lognormal'", numericInput("shift", "Shift", 200)),
      conditionalPanel("input.family == 'ordinal'", textInput("thresholds", "Thresholds", "-2, -0.6, 0.6, 2")),
      tags$hr(),
      actionButton("simulate", "Simulate", class = "btn-primary"),
      downloadButton("dl_spec", "Spec (.json)"), downloadButton("dl_data", "Data (.csv)")
    ),
    mainPanel(
      width = 8,
      tabsetPanel(
        tabPanel("Design spec (JSON)",
          p("Download this spec and run it unchanged in the R or Python package for identical data."),
          verbatimTextOutput("json")),
        tabPanel("Data", verbatimTextOutput("dims"), tableOutput("head")),
        tabPanel("Summary & plot", verbatimTextOutput("summary"), plotOutput("plot", height = "320px")),
        tabPanel("Power (Gaussian)",
          numericInput("n_sims", "Simulations", 500, min = 100, max = 2000, step = 100),
          actionButton("run_power", "Run"), verbatimTextOutput("power_out"))
      )
    )
  )
)

server <- function(input, output, session) {
  current_spec <- reactive(build_spec(list(
    name = input$name, seed = input$seed, n_subject = input$n_subject,
    include_items = input$include_items, n_item = input$n_item, design_kind = input$design_kind,
    factor_name = input$factor_name, lev1 = input$lev1, lev2 = input$lev2,
    intercept = input$intercept, effect = input$effect,
    subj_int_sd = input$subj_int_sd, subj_slope_sd = input$subj_slope_sd, subj_corr = input$subj_corr,
    item_int_sd = input$item_int_sd, item_slope_sd = input$item_slope_sd, item_corr = input$item_corr,
    family = input$family, resp_name = input$resp_name, sigma = input$sigma,
    shift = input$shift, thresholds = input$thresholds)))

  data <- eventReactive(input$simulate, simulate_design(current_spec()), ignoreNULL = FALSE)
  output$json <- renderText(spec_json(current_spec()))
  output$dims <- renderText({ d <- data(); sprintf("Simulated %d rows x %d columns.", nrow(d), ncol(d)) })
  output$head <- renderTable(head(data(), 10))
  output$summary <- renderPrint({
    d <- data(); yn <- current_spec()$response$name; fn <- input$factor_name
    if (is.numeric(d[[yn]])) print(do.call(data.frame, aggregate(d[[yn]], list(d[[fn]]),
        function(x) c(mean = mean(x), sd = sd(x), n = length(x))))) else print(table(d[[fn]], d[[yn]]))
  })
  output$plot <- renderPlot({
    d <- data(); yn <- current_spec()$response$name; fn <- input$factor_name
    if (is.numeric(d[[yn]])) boxplot(d[[yn]] ~ d[[fn]], xlab = fn, ylab = yn, col = c("#2C6FB0", "#B0402C"))
    else barplot(table(d[[fn]], d[[yn]]), beside = TRUE, legend = TRUE, col = c("#2C6FB0", "#B0402C"))
  })
  power_result <- reactiveVal(NULL)
  observeEvent(input$run_power, {
    spec <- current_spec()
    if (spec$response$family != "gaussian" || identical(input$design_kind, "within")) {
      power_result(list(msg = "This in-browser backend covers the two-group Gaussian design.\nInstall the package for crossed mixed-effects power.")); return()
    }
    power_result(power_design(spec, n_sims = min(as.integer(input$n_sims), 2000L)))
  })
  output$power_out <- renderPrint({
    r <- power_result(); if (is.null(r)) { cat("Click Run."); return() }
    if (!is.null(r$msg)) { cat(r$msg); return() }
    cat(sprintf("Power: %.3f | Type S: %.4f | Type M: %.3f | true %.3f, mean est %.3f\n",
                r$power, r$type_s, r$type_m, r$true_effect, r$mean_estimate))
  })
  output$dl_spec <- downloadHandler(function() paste0(input$name, ".json"),
    content = function(file) writeLines(spec_json(current_spec()), file))
  output$dl_data <- downloadHandler(function() paste0(input$name, ".csv"),
    content = function(file) write.csv(simulate_design(current_spec()), file, row.names = FALSE))
}

shinyApp(ui, server)
