# pilotr (lite), the serverless (shinylive / webR) build of the light path.
#
# Runs entirely in the browser via WebAssembly, with no server and no limit on the number
# of concurrent users. It covers design building (point-and-click, or pasting a spec),
# simulation, summaries and plots, two-group Gaussian simulation-based power (a point
# estimate and a power curve), a reproducible R script, and spec/data download. Crossed
# mixed-effects power (lme4) is provided by the installed pilotr package rather than this
# browser build. The engine and spec-builder files are staged next to this app by
# build_shinylive.R.

library(shiny)
library(bslib)
library(ggplot2)
for (f in c("core.R", "simulate.R", "parallel.R", "power.R", "spec_builder.R")) source(f)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (is.character(a) && !nzchar(a))) b else a

DOCS <- "https://pablobernabeu.github.io/pilotr/"
GH   <- "https://github.com/pablobernabeu/pilotr"

# Favicon embedded as a data URI so the serverless build stays self-contained (no separate asset).
FAVICON <- "data:image/png;base64,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"

# Reasonable intercept/effect (and noise) defaults per family, on each family's own scale.
# Switching family in the controls resets these so the point-and-click design stays valid
# (for example, a Gaussian-scale intercept of 100 would saturate a logit/log family).
FAMILY_DEFAULTS <- list(
  gaussian          = list(intercept = 100, effect = 5,   sigma = 10),
  shifted_lognormal = list(intercept = 6,   effect = 0.1, sigma = 0.3),
  bernoulli         = list(intercept = 0,   effect = 0.5),
  poisson           = list(intercept = 1.5, effect = 0.3),
  ordinal           = list(intercept = 0,   effect = 0.8),
  beta              = list(intercept = 0,   effect = 0.8, phi = 8)
)

PALETTE <- c("#2C6FB0", "#B0402C", "#2E8B57", "#8E6FB0", "#C8922A", "#5A5A5A")
pal <- function(k) if (k <= length(PALETTE)) PALETTE[seq_len(k)] else grDevices::hcl.colors(k, "Dynamic")

theme <- bs_theme(version = 5, primary = "#2C6FB0", "border-radius" = "0.5rem")

copy_btn <- function(target, label = "Copy")
  tags$button(label, class = "btn btn-sm btn-outline-secondary",
              onclick = sprintf("navigator.clipboard.writeText(document.getElementById('%s').innerText)", target))

# Immediate client-side feedback. webR runs the server in a worker, so the main thread can
# repaint on click while the computation is still in flight.
computing_js <- tags$script(HTML(
  "document.addEventListener('click', function(e){",
  "  var b = e.target.closest ? e.target.closest('#run_power, #run_curve') : null;",
  "  if (b) { var el = document.getElementById('power_out');",
  "           if (el) el.innerText = 'Computing… This runs in your browser and can take a few seconds.'; }",
  "}, true);"))

# Content-width, left-packed, striped data tables (to avoid the wide gap when a table has
# few columns), together with tighter sidebar rules.
app_css <- tags$style(HTML(
  "table.shiny-table { width: auto !important; align-self: flex-start; margin-bottom: 0; }",
  ".shiny-table th, .shiny-table td { padding: 0.25rem 0.9rem; }",
  ".sidebar hr { margin: 0.5rem 0; }",
  # Buttons: clear, smooth hover feedback (the solid primary did not visibly react).
  ".btn { transition: background-color .15s ease, border-color .15s ease, color .15s ease, box-shadow .15s ease; }",
  ".btn-primary:hover { background-color: #245d93; border-color: #21568a; }",
  ".btn-primary:active { background-color: #1f4f80 !important; border-color: #1f4f80 !important; }",
  ".btn-primary:hover, .btn-outline-primary:hover { box-shadow: 0 2px 7px rgba(44,111,176,.28); }",
  # Guide: section sub-headings and a tidy footer of links.
  ".guide-body h5 { margin-top: 1.6rem; margin-bottom: .55rem; font-size: 1.05rem; font-weight: 600;",
  "  color: #2C6FB0; border-bottom: 1px solid #e7eef5; padding-bottom: .3rem; }",
  ".guide-links { display: flex; flex-wrap: wrap; align-items: center; gap: .4rem .8rem;",
  "  margin-top: 1.6rem; padding-top: .9rem; border-top: 1px solid #e7eef5; font-size: .92rem; }",
  ".guide-links a { text-decoration: none; color: #2C6FB0; font-weight: 500; }",
  ".guide-links a:hover { text-decoration: underline; }",
  ".guide-links .sep { color: #c4ccd4; }"))

# ---------------------------------------------------------------------------------------------
guide <- card(
  card_header("How to use pilotr (lite)"),
  card_body(
    class = "guide-body",
    tags$p(tags$b("pilotr"), " simulates experimental and behavioural data from a portable ",
           "design specification. This browser version runs entirely on your device, with no ",
           "server and nothing uploaded. It covers the light path. You can build a design, ",
           "simulate data, inspect it, estimate two-group Gaussian power and export a ",
           "reproducible script or the specification itself."),
    tags$ol(
      tags$li("Describe the design in the left panel. This includes the sample sizes, the ",
              "factor and its two levels, the fixed intercept and effect, and a response ",
              "family. Changing the family resets the intercept and effect to reasonable ",
              "values for that family's scale. For within-subjects and crossed designs, you ",
              "can add random effects."),
      tags$li("Select ", tags$b("Simulate"), ", then read the tabs. These show the exact JSON ",
              "design specification, the simulated data, a summary and plot, a ",
              "simulation-based power analysis, and a reproducible R script."),
      tags$li("To take the design further, download or copy the specification (", tags$code(".json"),
              ") or the data (", tags$code(".csv"), "). You can also copy the R script and run ",
              "it in the installed package to reproduce the same data exactly.")
    ),
    tags$h5("Field guide"),
    tags$ul(
      tags$li(tags$b("Effect (±0.5)"), ", the difference between the two levels on ",
              "the response scale. This is the identity scale for Gaussian, the log scale ",
              "for RT and Poisson, and the logit scale for accuracy, ordinal and Beta."),
      tags$li(tags$b("Residual SD / Precision"), ", the Gaussian and shifted-lognormal noise ",
              "(an SD). Beta instead uses a precision (phi)."),
      tags$li(tags$b("Random-effect SDs"), ", the by-subject and by-item intercept and slope ",
              "standard deviations and their correlation, for crossed designs."),
      tags$li(tags$b("Seed"), ", which fixes the pseudo-random draw, so the same specification ",
              "and seed give identical data here, in R and in Python.")
    ),
    tags$h5("Advanced: paste a spec"),
    tags$p("The point-and-click controls cover common two-level designs across six response ",
           "families. The ", tags$b("Advanced: paste a full JSON spec"), " box in the sidebar ",
           "overrides those controls and gives access to the rest of the engine. This includes ",
           "continuous predictors and interactions, additional grouping factors (nesting), and ",
           "partial crossing. The example buttons there illustrate the format."),
    tags$div(
      class = "alert alert-info", role = "alert",
      tags$b("This is the lite build. "),
      "Crossed mixed-effects power analysis (via ", tags$code("lme4"),
      ") and precision/ROPE design analysis run in the installed ", tags$b("pilotr"),
      " R and Python packages rather than in the browser. The specification you build here ",
      "drives all three."
    ),
    tags$div(
      class = "guide-links",
      tags$a(href = DOCS, target = "_blank", "Documentation"),
      tags$span(class = "sep", "·"),
      tags$a(href = GH, target = "_blank", "Source and packages (R and Python)"),
      tags$span(class = "sep", "·"),
      tags$a(href = paste0(GH, "/blob/main/spec/SPEC.md"), target = "_blank", "Specification format")
    )
  )
)

ui <- page_sidebar(
  tags$head(tags$link(rel = "icon", type = "image/png", href = FAVICON)),
  title = "pilotr (lite): design and simulate in your browser",
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
      column(4, textInput("factor_name", "Factor", "group")),
      column(4, textInput("lev1", "Level 1", "control")),
      column(4, textInput("lev2", "Level 2", "treatment"))
    ),
    fluidRow(
      column(6, numericInput("intercept", "Intercept", 100)),
      column(6, numericInput("effect", "Effect (±0.5)", 5))
    ),
    selectInput("family", "Response family",
                c("Gaussian" = "gaussian", "Shifted lognormal (RT)" = "shifted_lognormal",
                  "Bernoulli (logit)" = "bernoulli", "Poisson (log)" = "poisson",
                  "Ordinal (Likert)" = "ordinal", "Beta (proportion)" = "beta")),
    textInput("resp_name", "Response name (blank = auto)", ""),
    conditionalPanel("input.family == 'gaussian' || input.family == 'shifted_lognormal'",
                     numericInput("sigma", "Residual SD", 10, min = 0)),
    conditionalPanel("input.family == 'shifted_lognormal'", numericInput("shift", "Shift", 200)),
    conditionalPanel("input.family == 'ordinal'",
                     textInput("thresholds", "Thresholds", "-2, -0.6, 0.6, 2")),
    conditionalPanel("input.family == 'beta'", numericInput("phi", "Precision (phi)", 8, min = 0.1)),
    tags$hr(),
    checkboxInput("use_pasted", tags$b("Advanced: paste a full JSON spec"), FALSE),
    conditionalPanel(
      "input.use_pasted == true",
      tags$small(class = "text-muted", "Overrides the controls above and gives access to the full engine."),
      textAreaInput("pasted", "Design spec (JSON)", "", height = "150px",
                    placeholder = "Paste a pilotr design spec (JSON)..."),
      tags$small(class = "text-muted",
                 "Need a starting point? Ready-to-run example specs are in ",
                 tags$a(href = paste0(GH, "/tree/main/spec/examples"), target = "_blank", "spec/examples"),
                 ", and the format is documented in ",
                 tags$a(href = paste0(GH, "/blob/main/spec/SPEC.md"), target = "_blank", "SPEC.md"), ".")
    ),
    tags$hr(),
    actionButton("simulate", "Simulate", class = "btn-primary w-100"),
    div(class = "mt-2 d-flex gap-2",
        downloadButton("dl_spec", "Spec (.json)", class = "btn-sm btn-outline-secondary"),
        downloadButton("dl_data", "Data (.csv)", class = "btn-sm btn-outline-secondary"))
  ),
  computing_js,
  app_css,
  navset_card_tab(
    id = "tabs",
    nav_panel("Guide", guide),
    nav_panel("Design spec",
      card(card_body(
        div(class = "d-flex justify-content-between align-items-center mb-2",
            tags$small("The portable source of truth. Run it unchanged in R or Python for identical data."),
            copy_btn("json", "Copy JSON")),
        verbatimTextOutput("json")))),
    nav_panel("Data",
      card(card_body(textOutput("dims"), tableOutput("head")))),
    nav_panel("Summary & plot",
      card(card_body(verbatimTextOutput("summary"), plotOutput("plot", height = "340px")))),
    nav_panel("Power",
      card(card_body(
        fluidRow(
          column(4, numericInput("n_sims", "Simulations", 300, min = 100, max = 3000, step = 100)),
          column(8, div(class = "mt-4 d-flex gap-2 flex-wrap",
                        actionButton("run_power", "Estimate power", class = "btn-primary"),
                        actionButton("run_curve", "Power curve over N", class = "btn-outline-primary")))
        ),
        tags$small(class = "text-muted", "Two-group Gaussian, in-browser. The curve uses up to 200 sims/point to stay responsive."),
        verbatimTextOutput("power_out"),
        plotOutput("power_plot", height = "320px")))),
    nav_panel("R script",
      card(card_body(
        div(class = "d-flex justify-content-between align-items-center mb-2",
            tags$small("A self-contained script that reproduces this design with the installed pilotr R package."),
            div(class = "d-flex gap-2",
                copy_btn("rscript", "Copy"),
                downloadButton("dl_script", "Download .R", class = "btn-sm btn-outline-secondary"))),
        verbatimTextOutput("rscript"))))
  )
)

# ---------------------------------------------------------------------------------------------
server <- function(input, output, session) {

  # Reset intercept/effect (and noise) to family-appropriate values when the family changes.
  observeEvent(input$family, {
    d <- FAMILY_DEFAULTS[[input$family]]
    if (is.null(d)) return()
    updateNumericInput(session, "intercept", value = d$intercept)
    updateNumericInput(session, "effect",    value = d$effect)
    if (!is.null(d$sigma)) updateNumericInput(session, "sigma", value = d$sigma)
    if (!is.null(d$phi))   updateNumericInput(session, "phi",   value = d$phi)
  }, ignoreInit = TRUE)

  parse_error <- reactiveVal(NULL)

  # Returns the spec, or NULL if a pasted spec is invalid (with the message in parse_error()).
  current_spec <- reactive({
    if (isTRUE(input$use_pasted) && nzchar(trimws(input$pasted %||% ""))) {
      spec <- tryCatch({
        tf <- tempfile(fileext = ".json"); writeLines(input$pasted, tf); load_spec(tf)
      }, error = function(e) e)
      if (inherits(spec, "error")) { parse_error(conditionMessage(spec)); return(NULL) }
      parse_error(NULL); spec
    } else {
      parse_error(NULL)
      build_spec(list(
        name = input$name, seed = input$seed, n_subject = input$n_subject,
        include_items = input$include_items, n_item = input$n_item, design_kind = input$design_kind,
        factor_name = input$factor_name, lev1 = input$lev1, lev2 = input$lev2,
        intercept = input$intercept, effect = input$effect,
        subj_int_sd = input$subj_int_sd, subj_slope_sd = input$subj_slope_sd, subj_corr = input$subj_corr,
        item_int_sd = input$item_int_sd, item_slope_sd = input$item_slope_sd, item_corr = input$item_corr,
        family = input$family, resp_name = input$resp_name, sigma = input$sigma,
        shift = input$shift, thresholds = input$thresholds, phi = input$phi))
    }
  })

  spec_req <- function() {
    s <- current_spec()
    validate(need(!is.null(s), parse_error() %||% "Please enter a valid design specification."))
    s
  }
  nsims <- function() {
    v <- suppressWarnings(as.integer(input$n_sims))
    if (length(v) == 0 || is.na(v)) 500L else max(100L, min(v, 3000L))
  }
  resp_name  <- function(spec) spec$response$name %||% "y"
  group_name <- function(spec) {
    f <- spec$factors
    if (length(f) >= 1 && !is.null(f[[1]]$name)) f[[1]]$name else NA_character_
  }

  sim_data <- eventReactive(input$simulate, {
    s <- current_spec(); if (is.null(s)) return(NULL); simulate_design(s)
  }, ignoreNULL = FALSE)

  data_req <- function() {
    d <- sim_data()
    validate(need(!is.null(d), parse_error() %||% "Please correct the specification, then select Simulate."))
    d
  }

  output$json <- renderText(spec_json(spec_req()))

  output$dims <- renderText({
    d <- data_req(); sprintf("Simulated %d rows x %d columns.", nrow(d), ncol(d))
  })
  output$head <- renderTable(head(data_req(), 12), striped = TRUE, hover = TRUE, spacing = "xs")

  output$summary <- renderPrint({
    d <- data_req(); yn <- resp_name(current_spec()); gn <- group_name(current_spec())
    y <- d[[yn]]
    if (is.na(gn) || is.null(d[[gn]])) {
      if (is.numeric(y)) print(round(c(mean = mean(y), sd = sd(y), min = min(y), max = max(y)), 3))
      else print(table(y))
    } else if (is.numeric(y)) {
      agg <- aggregate(y, setNames(list(d[[gn]]), gn),
                       function(x) round(c(mean = mean(x), sd = sd(x), n = length(x)), 3))
      print(data.frame(agg[gn], agg$x, check.names = FALSE), row.names = FALSE)
    } else {
      print(table(d[[gn]], y))
    }
  })

  output$plot <- renderPlot({
    d <- data_req(); yn <- resp_name(current_spec()); gn <- group_name(current_spec())
    y <- d[[yn]]; base <- theme_minimal(base_size = 14)
    has_grp <- !is.na(gn) && !is.null(d[[gn]])
    if (has_grp && is.numeric(y)) {
      ggplot(d, aes(.data[[gn]], .data[[yn]], fill = .data[[gn]])) +
        geom_boxplot(alpha = 0.85, outlier.alpha = 0.35) +
        scale_fill_manual(values = pal(2), guide = "none") +
        labs(x = gn, y = yn, title = paste(yn, "by", gn)) + base
    } else if (has_grp) {
      ggplot(d, aes(.data[[yn]], fill = .data[[gn]])) +
        geom_bar(position = "dodge") +
        scale_fill_manual(values = pal(length(unique(d[[gn]]))), name = gn) +
        labs(x = yn, y = "count", title = paste(yn, "by", gn)) + base
    } else if (is.numeric(y)) {
      ggplot(d, aes(.data[[yn]])) +
        geom_histogram(bins = 30, fill = PALETTE[1], colour = "white") +
        labs(x = yn, y = "count", title = paste("Distribution of", yn)) + base
    } else {
      ggplot(d, aes(.data[[yn]])) +
        geom_bar(fill = PALETTE[1]) +
        labs(x = yn, y = "count", title = paste("Distribution of", yn)) + base
    }
  })

  # ---- power (two-group Gaussian) ----
  gaussian_two_group <- function(spec) {
    identical(spec$response$family, "gaussian") &&
      length(spec$factors) >= 1 && !is.null(spec$factors[[1]]$between)
  }
  not_supported_msg <- paste0(
    "The in-browser demo runs power only for the two-group Gaussian design.\n",
    "This design needs the installed package. Download the spec with the Spec (.json)\n",
    "button, then run one of the following.\n\n",
    "R (crossed mixed-effects power and a power curve, via lme4):\n",
    "    library(pilotr)\n",
    "    spec <- load_spec(\"design.json\")\n",
    "    power_mixed(spec, n_sims = 200)\n",
    "    power_curve_mixed(spec, subject_ns = c(20, 40, 60), n_sims = 200)\n\n",
    "Python (statsmodels backend):\n",
    "    from pilotr import load_spec, power_mixed\n",
    "    power_mixed(load_spec(\"design.json\"), n_sims=200)\n\n",
    "Install:\n",
    "    remotes::install_github(\"pablobernabeu/pilotr\", subdir = \"r/pilotr\")   # R\n",
    "    pip install \"git+https://github.com/pablobernabeu/pilotr#subdirectory=python\"   # Python")

  power_out  <- reactiveVal("Select “Estimate power” to run a simulation-based power analysis.")
  power_plot <- reactiveVal(NULL)

  observeEvent(input$run_power, {
    spec <- current_spec()
    if (is.null(spec)) { power_out(parse_error() %||% "Please correct the specification first."); return() }
    if (!gaussian_two_group(spec)) { power_out(not_supported_msg); power_plot(NULL); return() }
    ns <- nsims(); r <- power_design(spec, n_sims = ns)
    power_out(sprintf(
      "Power: %.3f   |   Type S: %.4f   |   Type M: %.3f\nTrue effect: %.3f   |   mean estimate: %.3f   (n_sims = %d)",
      r$power, r$type_s, r$type_m, r$true_effect, r$mean_estimate, ns))
  })

  observeEvent(input$run_curve, {
    spec <- current_spec()
    if (is.null(spec)) { power_out(parse_error() %||% "Please correct the specification first."); power_plot(NULL); return() }
    if (!gaussian_two_group(spec)) { power_out(not_supported_msg); power_plot(NULL); return() }
    base_n <- spec$units$subject$n
    grid <- unique(round(base_n * c(0.5, 0.75, 1, 1.5, 2)))
    grid <- grid[grid >= 4]
    ns   <- min(nsims(), 200L)   # keep the in-browser sweep responsive
    pw   <- vapply(grid, function(n) {
      s <- spec; s$units$subject$n <- as.integer(n); power_design(s, n_sims = ns)$power
    }, numeric(1))
    power_out(sprintf("Power curve at n_sims = %d per point. N subjects = %s.",
                      ns, paste(grid, collapse = ", ")))
    power_plot(list(grid = grid, pw = pw))
  })

  # The script above shows "Computing…" the instant a power button is clicked. Clearing it
  # needs the server to push a fresh value, but Shiny skips sending an output whose value is
  # unchanged, so a repeated (identical) result would otherwise leave "Computing…" stuck.
  # Toggling an invisible trailing space by the click count makes each value distinct, so the
  # update is always sent. The space does not show in the monospaced output.
  output$power_out <- renderText({
    n <- sum(input$run_power, input$run_curve)
    paste0(power_out(), strrep(" ", n %% 2))
  })
  output$power_plot <- renderPlot({
    pc <- power_plot(); if (is.null(pc)) return(NULL)
    df <- data.frame(n = pc$grid, power = pc$pw)
    ggplot(df, aes(n, power)) +
      geom_hline(yintercept = 0.8, linetype = 2, colour = "#888888") +
      annotate("text", x = min(df$n), y = 0.8, label = "0.80 target",
               hjust = 0, vjust = -0.6, colour = "#888888", size = 3.6) +
      geom_line(colour = PALETTE[1], linewidth = 0.9) +
      geom_point(colour = PALETTE[1], size = 3) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(x = expression(italic(N) ~ "subjects"), y = "Power", title = "Power curve") +
      theme_minimal(base_size = 14)
  })

  # ---- reproducible R script ----
  output$rscript <- renderText(generate_r_script(spec_req()))
  output$dl_script <- downloadHandler(
    function() paste0(current_spec()$name %||% "design", ".R"),
    content = function(file) {
      s <- current_spec()
      writeLines(if (is.null(s)) "# Invalid specification. Please correct the pasted JSON." else generate_r_script(s), file)
    })

  # ---- downloads (degrade gracefully on an invalid pasted spec) ----
  output$dl_spec <- downloadHandler(
    function() paste0(current_spec()$name %||% "design", ".json"),
    content = function(file) {
      s <- current_spec()
      writeLines(if (is.null(s)) "{}" else spec_json(s), file)
    })
  output$dl_data <- downloadHandler(
    function() paste0(current_spec()$name %||% "design", ".csv"),
    content = function(file) {
      s <- current_spec()
      if (is.null(s)) writeLines("", file) else write.csv(simulate_design(s), file, row.names = FALSE)
    })

  # Keep the download links live even while their tab is hidden, so a click always works
  # (Shiny suspends hidden outputs by default, which left the R-script download href empty).
  for (id in c("dl_script", "dl_spec", "dl_data"))
    outputOptions(output, id, suspendWhenHidden = FALSE)
}

shinyApp(ui, server)
