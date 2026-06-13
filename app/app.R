# simdgp no-code app -- the third interface over the shared design spec.
#
# This is a THIN client: it builds the portable JSON spec from point-and-click inputs and
# calls the simdgp R package to simulate and to run power analysis. The same downloaded
# spec runs identically in the R and Python packages -- "one model, three interfaces".
#
# Run with:  shiny::runApp("toolkit/app")

library(shiny)
library(jsonlite)

# ---- locate and source the simdgp package sources + the spec builder ----
.find <- function(rel) {
  for (base in c(".", "..", "../..", "toolkit", normalizePath(".", mustWork = FALSE))) {
    p <- file.path(base, rel)
    if (file.exists(p)) return(normalizePath(p))
  }
  stop("could not locate ", rel)
}
source(.find("app/spec_builder.R"))
for (f in c("core.R", "simulate.R", "power.R", "power_mixed.R"))
  source(.find(file.path("r", "simdgp", "R", f)))

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
        column(6, numericInput("intercept", "Intercept (grand mean / log-rate / logit)", 100)),
        column(6, numericInput("effect", "Effect size (coef. on -0.5/+0.5 contrast)", 5))
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
        tabPanel("Data",
          verbatimTextOutput("dims"), tableOutput("head")),
        tabPanel("Summary & plot",
          verbatimTextOutput("summary"), plotOutput("plot", height = "320px")),
        tabPanel("Power & design analysis",
          p("Simulation-based power with Type S / Type M (Gelman & Carlin, 2014)."),
          numericInput("n_sims", "Simulations", 1000, min = 100, step = 100),
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

  data <- eventReactive(input$simulate, { simulate_design(current_spec()) }, ignoreNULL = FALSE)

  output$json <- renderText(spec_json(current_spec()))

  output$dims <- renderText({
    d <- data(); sprintf("Simulated %d rows x %d columns (seed %d).", nrow(d), ncol(d), input$seed)
  })
  output$head <- renderTable(head(data(), 10))

  output$summary <- renderPrint({
    d <- data(); yn <- current_spec()$response$name; fn <- input$factor_name
    if (is.numeric(d[[yn]])) {
      agg <- aggregate(d[[yn]], list(d[[fn]]), function(x) c(mean = mean(x), sd = sd(x), n = length(x)))
      cat("Mean (SD) of", yn, "by", fn, ":\n"); print(do.call(data.frame, agg))
    } else {
      cat("Counts of", yn, "by", fn, ":\n"); print(table(d[[fn]], d[[yn]]))
    }
  })

  output$plot <- renderPlot({
    d <- data(); yn <- current_spec()$response$name; fn <- input$factor_name
    if (is.numeric(d[[yn]])) {
      boxplot(d[[yn]] ~ d[[fn]], xlab = fn, ylab = yn, col = c("#2C6FB0", "#B0402C"),
              main = paste("Distribution of", yn))
    } else {
      barplot(table(d[[fn]], d[[yn]]), beside = TRUE, legend = TRUE,
              col = c("#2C6FB0", "#B0402C"), main = paste("Counts of", yn))
    }
  })

  output$repro <- renderText({
    sprintf(paste0(
      "# 1. Download the spec from the first tab as design.json\n\n",
      "# --- R ---\n",
      "library(simdgp)\n",
      "d <- simulate_design(\"design.json\")\n\n",
      "# --- Python ---\n",
      "from simdgp import simulate\n",
      "d = simulate(\"design.json\")\n",
      "d.to_csv(\"data.csv\")\n\n",
      "# Same spec + seed => identical data in all three interfaces."))
  })

  observeEvent(input$run_power, {
    spec <- current_spec()
    output$power_out <- renderPrint({
      if (spec$response$family != "gaussian" || identical(input$design_kind, "within")) {
        cat("The in-app power backend covers the two-group Gaussian design.\n")
        cat("For crossed mixed-effects designs use power_mixed() in the R package\n",
            "(fits lme4; can take minutes) or power_curve_mixed() for a curve.\n")
        return(invisible())
      }
      res <- power_design(spec, n_sims = input$n_sims)
      cat(sprintf("Simulations : %d\n", res$n_sims))
      cat(sprintf("Power       : %.3f\n", res$power))
      cat(sprintf("Type S error: %.4f\n", res$type_s))
      cat(sprintf("Type M (exag): %.3f\n", res$type_m))
      cat(sprintf("True effect  : %.3f  | mean estimate: %.3f\n", res$true_effect, res$mean_estimate))
    })
  })

  output$dl_spec <- downloadHandler(
    filename = function() paste0(input$name, ".json"),
    content = function(file) writeLines(spec_json(current_spec()), file))
  output$dl_data <- downloadHandler(
    filename = function() paste0(input$name, "_seed", input$seed, ".csv"),
    content = function(file) write.csv(simulate_design(current_spec()), file, row.names = FALSE))
}

shinyApp(ui, server)
