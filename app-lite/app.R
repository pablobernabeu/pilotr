# pilotr LITE -- the serverless (shinylive / webR) build of the light path.
#
# Runs entirely in the browser via WebAssembly: zero server, unlimited concurrent users.
# Covers design building (point-and-click OR paste any spec), simulation, summaries and
# plots, two-group Gaussian simulation-based power (point estimate + a power curve), a
# reproducible R script, and spec/data download. Heavy crossed mixed-effects power (lme4)
# is intentionally a feature of the installed pilotr package, not this browser build.
# The engine + spec-builder files are staged next to this app by build_shinylive.R.

library(shiny)
library(bslib)
for (f in c("core.R", "simulate.R", "power.R", "spec_builder.R")) source(f)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (is.character(a) && !nzchar(a))) b else a

DOCS   <- "https://pablobernabeu.github.io/pilotr/"
GH     <- "https://github.com/pablobernabeu/pilotr"

# ---- worked examples for the "paste a spec" advanced mode (the full engine) ----------------
EX_BETA <- '{
  "name": "beta_proportion_demo",
  "seed": 2024,
  "units": { "subject": { "n": 120 } },
  "factors": [
    { "name": "group", "levels": ["control", "treatment"],
      "contrasts": { "effect": [-0.5, 0.5] }, "between": "subject" }
  ],
  "fixed": { "intercept": 0.0, "coefficients": { "effect": 0.8 } },
  "response": { "family": "beta", "name": "proportion", "phi": 8 }
}'

EX_CONTINUOUS <- '{
  "name": "continuous_predictor_demo",
  "seed": 2024,
  "units": { "subject": { "n": 40 }, "item": { "n": 30 } },
  "predictors": [ { "name": "frequency", "varies_by": "item", "mean": 0, "sd": 1 } ],
  "factors": [],
  "fixed": { "intercept": 6.0, "coefficients": { "frequency": -0.05 } },
  "random": {
    "subject": { "intercept_sd": 0.12, "slopes": { "frequency": 0.04 },
                 "correlations": { "intercept,frequency": 0.2 } },
    "item": { "intercept_sd": 0.08 }
  },
  "response": { "family": "shifted_lognormal", "name": "RT", "sigma": 0.3, "shift": 200, "round": 4 }
}'

PALETTE <- c("#2C6FB0", "#B0402C", "#2E8B57", "#8E6FB0")

theme <- bs_theme(version = 5, primary = "#2C6FB0", "border-radius" = "0.5rem")

# ---------------------------------------------------------------------------------------------
guide <- card(
  card_header("How to use pilotr (lite)"),
  card_body(
    tags$p(tags$b("pilotr"), " simulates experimental and behavioral data from a portable ",
           "design specification. This browser version runs entirely on your device (no ",
           "server, nothing uploaded) and covers the light path: build a design, simulate ",
           "data, inspect it, estimate two-group Gaussian power, and export a reproducible ",
           "script or the spec itself."),
    tags$ol(
      tags$li(tags$b("Describe the design"), " in the left panel — sample sizes, the factor ",
              "and its two levels, the fixed intercept and effect, and a response family. ",
              "For within/crossed designs you can add by-subject and by-item random effects."),
      tags$li(tags$b("Click Simulate"), ", then read the tabs: the exact JSON ",
              tags$b("Design spec"), ", the simulated ", tags$b("Data"), ", a ",
              tags$b("Summary & plot"), ", simulation-based ", tags$b("Power"), ", and a ",
              tags$b("Reproducible R script"), "."),
      tags$li(tags$b("Take it further."), " Download the spec (", tags$code(".json"),
              ") or data (", tags$code(".csv"), "), or copy the R script and run it in the ",
              "installed package to reproduce the same data bit-for-bit.")
    ),
    tags$h6("Field guide"),
    tags$ul(
      tags$li(tags$b("Effect (on -0.5/+0.5)"), " — the difference between the two levels on ",
              "the response scale (log scale for RT/Poisson, logit for accuracy/ordinal)."),
      tags$li(tags$b("Residual SD"), " — the within-cell noise (Gaussian / shifted-lognormal)."),
      tags$li(tags$b("Random-effect SDs"), " — by-subject / by-item intercept and slope ",
              "standard deviations, and their correlation, for crossed designs."),
      tags$li(tags$b("Seed"), " — fixes the pseudo-random draw, so the same spec + seed gives ",
              "identical data here, in R, and in Python.")
    ),
    tags$h6("Advanced: paste a spec"),
    tags$p("The point-and-click controls cover common two-level designs. The ",
           tags$b("Advanced: paste a full JSON spec"), " box in the sidebar unlocks the full ",
           "engine — Beta proportions, continuous predictors and interactions, additional ",
           "grouping factors (nesting), and partial crossing. Use the example buttons there ",
           "to see the format."),
    tags$div(
      class = "alert alert-info", role = "alert",
      tags$b("This is the lite build. "),
      "Crossed mixed-effects power analysis (via ", tags$code("lme4"),
      ") and precision/ROPE design analysis run in the installed ", tags$b("pilotr"),
      " R and Python packages, not in the browser. The spec you build here drives all three."
    ),
    tags$p(
      tags$a(href = DOCS, target = "_blank", "Documentation"), " · ",
      tags$a(href = GH, target = "_blank", "Source & packages (R / Python)"), " · ",
      tags$a(href = paste0(GH, "/blob/main/spec/SPEC.md"), target = "_blank", "Specification format")
    )
  )
)

ui <- page_sidebar(
  title = "pilotr (lite) — design & simulate in your browser",
  theme = theme,
  sidebar = sidebar(
    width = 360,
    title = "Design",
    textInput("name", "Design name", "my_design"),
    numericInput("seed", "Seed", 2024, step = 1),
    fluidRow(
      column(7, selectInput("design_kind", "Design",
                            c("Between-subjects" = "between", "Within / crossed" = "within"))),
      column(5, numericInput("n_subject", "N subjects", 64, min = 2))
    ),
    conditionalPanel(
      "input.design_kind == 'within'",
      checkboxInput("include_items", "Crossed with items", TRUE),
      conditionalPanel("input.include_items", numericInput("n_item", "N items", 24, min = 2)),
      tags$small(tags$b("By-subject random effects")),
      fluidRow(
        column(4, numericInput("subj_int_sd", "Int SD", 0.12, min = 0)),
        column(4, numericInput("subj_slope_sd", "Slope SD", 0.04, min = 0)),
        column(4, numericInput("subj_corr", "corr", 0.2, min = -1, max = 1))
      ),
      conditionalPanel(
        "input.include_items",
        tags$small(tags$b("By-item random effects")),
        fluidRow(
          column(4, numericInput("item_int_sd", "Int SD", 0.08, min = 0)),
          column(4, numericInput("item_slope_sd", "Slope SD", 0.02, min = 0)),
          column(4, numericInput("item_corr", "corr", -0.1, min = -1, max = 1))
        )
      )
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
    conditionalPanel("input.family == 'ordinal'",
                     textInput("thresholds", "Thresholds", "-2, -0.6, 0.6, 2")),
    tags$hr(),
    checkboxInput("use_pasted", tags$b("Advanced: paste a full JSON spec"), FALSE),
    conditionalPanel(
      "input.use_pasted == true",
      tags$small("Overrides the controls above. Unlocks the full engine."),
      textAreaInput("pasted", NULL, "", height = "160px",
                    placeholder = "Paste a pilotr design spec (JSON)..."),
      tags$small("Load example: "),
      actionLink("ex_beta", "Beta"), " · ",
      actionLink("ex_cont", "continuous predictor")
    ),
    tags$hr(),
    actionButton("simulate", "Simulate", class = "btn-primary w-100"),
    div(class = "mt-2 d-flex gap-2",
        downloadButton("dl_spec", "Spec (.json)", class = "btn-sm btn-outline-secondary"),
        downloadButton("dl_data", "Data (.csv)", class = "btn-sm btn-outline-secondary"))
  ),
  navset_card_tab(
    id = "tabs",
    nav_panel("Guide", guide),
    nav_panel("Design spec",
      card(card_body(
        tags$p("This portable spec is the source of truth. Download it and run it unchanged ",
               "in the R or Python package for identical data."),
        verbatimTextOutput("json")))),
    nav_panel("Data",
      card(card_body(textOutput("dims"), tableOutput("head")))),
    nav_panel("Summary & plot",
      card(card_body(verbatimTextOutput("summary"), plotOutput("plot", height = "340px")))),
    nav_panel("Power",
      card(card_body(
        fluidRow(
          column(4, numericInput("n_sims", "Simulations", 500, min = 100, max = 3000, step = 100)),
          column(8, div(class = "mt-4",
                        actionButton("run_power", "Estimate power", class = "btn-primary"),
                        actionButton("run_curve", "Power curve over N", class = "btn-outline-primary")))
        ),
        verbatimTextOutput("power_out"),
        plotOutput("power_plot", height = "320px")))),
    nav_panel("R script",
      card(card_body(
        tags$p("A self-contained script that reproduces this design with the installed ",
               tags$b("pilotr"), " R package."),
        downloadButton("dl_script", "Download .R", class = "btn-sm btn-outline-secondary mb-2"),
        verbatimTextOutput("rscript"))))
  )
)

# ---------------------------------------------------------------------------------------------
server <- function(input, output, session) {

  observeEvent(input$ex_beta, {
    updateCheckboxInput(session, "use_pasted", value = TRUE)
    updateTextAreaInput(session, "pasted", value = EX_BETA)
  })
  observeEvent(input$ex_cont, {
    updateCheckboxInput(session, "use_pasted", value = TRUE)
    updateTextAreaInput(session, "pasted", value = EX_CONTINUOUS)
  })

  current_spec <- reactive({
    if (isTRUE(input$use_pasted) && nzchar(trimws(input$pasted %||% ""))) {
      spec <- tryCatch({
        tf <- tempfile(fileext = ".json"); writeLines(input$pasted, tf); load_spec(tf)
      }, error = function(e) e)
      validate(need(!inherits(spec, "error"),
                    paste("Could not parse the spec:", conditionMessage(spec))))
      spec
    } else {
      build_spec(list(
        name = input$name, seed = input$seed, n_subject = input$n_subject,
        include_items = input$include_items, n_item = input$n_item, design_kind = input$design_kind,
        factor_name = input$factor_name, lev1 = input$lev1, lev2 = input$lev2,
        intercept = input$intercept, effect = input$effect,
        subj_int_sd = input$subj_int_sd, subj_slope_sd = input$subj_slope_sd, subj_corr = input$subj_corr,
        item_int_sd = input$item_int_sd, item_slope_sd = input$item_slope_sd, item_corr = input$item_corr,
        family = input$family, resp_name = input$resp_name, sigma = input$sigma,
        shift = input$shift, thresholds = input$thresholds))
    }
  })

  resp_name  <- function(spec) spec$response$name %||% "y"
  group_name <- function(spec) {
    f <- spec$factors
    if (length(f) >= 1 && !is.null(f[[1]]$name)) f[[1]]$name else NA_character_
  }

  sim_data <- eventReactive(input$simulate, simulate_design(current_spec()), ignoreNULL = FALSE)

  output$json <- renderText(spec_json(current_spec()))

  output$dims <- renderText({
    d <- sim_data(); sprintf("Simulated %d rows x %d columns.", nrow(d), ncol(d))
  })
  output$head <- renderTable(head(sim_data(), 12))

  output$summary <- renderPrint({
    d <- sim_data(); yn <- resp_name(current_spec()); gn <- group_name(current_spec())
    y <- d[[yn]]
    if (is.na(gn) || is.null(d[[gn]])) {
      if (is.numeric(y)) print(round(c(mean = mean(y), sd = sd(y), min = min(y), max = max(y)), 3))
      else print(table(y))
    } else if (is.numeric(y)) {
      print(do.call(data.frame, aggregate(y, list(d[[gn]]),
            function(x) round(c(mean = mean(x), sd = sd(x), n = length(x)), 3))))
    } else {
      print(table(d[[gn]], y))
    }
  })

  output$plot <- renderPlot({
    d <- sim_data(); yn <- resp_name(current_spec()); gn <- group_name(current_spec())
    y <- d[[yn]]; op <- par(mar = c(4, 4, 1, 1)); on.exit(par(op))
    if (!is.na(gn) && !is.null(d[[gn]]) && is.numeric(y)) {
      boxplot(y ~ d[[gn]], xlab = gn, ylab = yn, col = PALETTE[1:2], border = "#333333")
    } else if (!is.na(gn) && !is.null(d[[gn]])) {
      barplot(table(d[[gn]], y), beside = TRUE, legend = TRUE, xlab = yn, col = PALETTE[1:2])
    } else if (is.numeric(y)) {
      hist(y, breaks = 30, col = PALETTE[1], border = "white", main = NULL, xlab = yn)
    } else {
      barplot(table(y), col = PALETTE[1], ylab = "count", xlab = yn)
    }
  })

  # ---- power (two-group Gaussian) ----
  gaussian_two_group <- function(spec) {
    identical(spec$response$family, "gaussian") &&
      length(spec$factors) >= 1 && !is.null(spec$factors[[1]]$between)
  }
  not_supported_msg <- paste0(
    "The in-browser backend estimates power for the two-group Gaussian design.\n",
    "For crossed mixed-effects power (lme4) and precision/ROPE analysis, install the\n",
    "pilotr package and run, e.g.:\n\n",
    "    library(pilotr)\n",
    "    pow <- power_mixed(load_spec(\"design.json\"), n_sims = 200)")

  power_out  <- reactiveVal("Click “Estimate power” to run a simulation-based power analysis.")
  power_plot <- reactiveVal(NULL)

  observeEvent(input$run_power, {
    spec <- current_spec()
    if (!gaussian_two_group(spec)) { power_out(not_supported_msg); return() }
    r <- power_design(spec, n_sims = min(as.integer(input$n_sims), 3000L))
    power_out(sprintf(
      "Power: %.3f   |   Type S: %.4f   |   Type M: %.3f\nTrue effect: %.3f   |   mean estimate: %.3f   (n_sims = %d)",
      r$power, r$type_s, r$type_m, r$true_effect, r$mean_estimate, min(as.integer(input$n_sims), 3000L)))
  })

  observeEvent(input$run_curve, {
    spec <- current_spec()
    if (!gaussian_two_group(spec)) { power_out(not_supported_msg); power_plot(NULL); return() }
    base_n <- spec$units$subject$n
    grid <- unique(round(base_n * c(0.5, 0.75, 1, 1.5, 2)))
    grid <- grid[grid >= 4]
    ns   <- min(as.integer(input$n_sims), 1500L)
    pw   <- vapply(grid, function(n) {
      s <- spec; s$units$subject$n <- as.integer(n); power_design(s, n_sims = ns)$power
    }, numeric(1))
    power_out(sprintf("Power curve at n_sims = %d per point. Total N = %s.",
                      ns, paste(grid, collapse = ", ")))
    power_plot(list(grid = grid, pw = pw))
  })

  output$power_out <- renderText(power_out())
  output$power_plot <- renderPlot({
    pc <- power_plot(); if (is.null(pc)) return(NULL)
    op <- par(mar = c(4, 4, 1, 1)); on.exit(par(op))
    plot(pc$grid, pc$pw, type = "b", pch = 19, col = PALETTE[1], lwd = 2,
         ylim = c(0, 1), xlab = "Total N", ylab = "Power")
    abline(h = 0.8, lty = 2, col = "#888888")
  })

  # ---- reproducible R script ----
  output$rscript <- renderText(generate_r_script(current_spec()))
  output$dl_script <- downloadHandler(
    function() paste0(current_spec()$name %||% "design", ".R"),
    content = function(file) writeLines(generate_r_script(current_spec()), file))

  # ---- downloads ----
  output$dl_spec <- downloadHandler(
    function() paste0(current_spec()$name %||% "design", ".json"),
    content = function(file) writeLines(spec_json(current_spec()), file))
  output$dl_data <- downloadHandler(
    function() paste0(current_spec()$name %||% "design", ".csv"),
    content = function(file) write.csv(simulate_design(current_spec()), file, row.names = FALSE))
}

shinyApp(ui, server)
