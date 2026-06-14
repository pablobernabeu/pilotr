# pilotr no-code app -- the third interface over the shared design spec.
#
# Thin client: every control writes into the portable JSON spec, which downloads and runs
# unchanged in the R and Python packages. Launch with pilotr::run_app() (installed) or
# shiny::runApp("r/pilotr/inst/app") (from source).

library(shiny)

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

# ---------------------------------------------------------------- UI ----
ui <- fluidPage(
  titlePanel("pilotr — design · simulate · power (one spec, three interfaces)"),
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
        tabPanel("Reproducible R script",
          p("Your no-code design as a self-contained R script. Download it, or ",
            tags$b("verify"), " it reproduces this exact data by running it in a clean R session."),
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
      shift = input$shift, thresholds = input$thresholds
    ))
  })

  data <- eventReactive(input$simulate, simulate_design(current_spec()), ignoreNULL = FALSE)

  output$json <- renderText(spec_json(current_spec()))
  output$dims <- renderText({ d <- data(); sprintf("Simulated %d rows x %d columns (seed %d).", nrow(d), ncol(d), input$seed) })
  output$head <- renderTable(head(data(), 10))

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
    if (is.numeric(d[[yn]])) {
      if (has_grp) boxplot(d[[yn]] ~ d[[fn]], xlab = fn, ylab = yn, col = c("#2C6FB0", "#B0402C"),
                           main = paste("Distribution of", yn))
      else hist(d[[yn]], col = "#2C6FB0", xlab = yn, main = paste("Distribution of", yn))
    } else if (has_grp) barplot(table(d[[fn]], d[[yn]]), beside = TRUE, legend = TRUE,
                                col = c("#2C6FB0", "#B0402C"), main = paste("Counts of", yn))
    else barplot(table(d[[yn]]), col = "#2C6FB0", main = paste("Counts of", yn))
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
    if (is.null(r)) return("Click 'Verify' to run the script in a fresh R process and confirm it reproduces this data.")
    if (!is.null(r$msg)) return(r$msg)
    if (isTRUE(r$ok)) sprintf("OK - reproduces identically in a clean R session.\n  %d rows (app) = %d rows (clean run); response checksum matches.", r$ref_n, r$n)
    else sprintf("MISMATCH - app %d rows vs clean run %d rows, or checksum differs.", r$ref_n, r$n)
  })

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
      p <- promises::future_promise({ pilotr::power_design(spec, n_sims = n) }, seed = TRUE)
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
