##
## This was a Shiny web application providing a simple interface to the
## simulation capabilities of the 'precautionary' package.
## Now I want to obtain an 'MCSE-free' version that generates
## safety schematics such as Fig 3 of Norris 2020c.
##
## TODO:
## 1. Refine event handling for custom doses
## 2. Run longer CPEs with pop-up progress bar?
## MARK III (?)
## 1. Obtain the safety schematic
## 2. Can I 'mark the spot' on schematic selected by hyperprior?
##    - Or is it perhaps a *region* on the schematic?

library(shiny)
library(shinyjs)
library(precautionary)
library(kableExtra)

ui <- fluidPage(
  shinyjs::useShinyjs(),
  shinyFeedback::useShinyFeedback(),

  singleton(tags$head(tags$link(rel="stylesheet", type = "text/css", href = "introjs.css"))),
  singleton(tags$head(tags$script(src="intro.js"))),
  tags$head(singleton(tags$script(src="events.js"))),

  singleton(tags$head(tags$link(rel="stylesheet", type = "text/css", href = "tweaks.css"))),

  ## Application title
  titlePanel("Predict Risks of High-Grade Toxicities in Dose-Escalation Trials"),

  ## Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      tags$fieldset(id="dose-levels",
      tags$legend("prespecified dose levels"),
      splitLayout(
        textInput(inputId = "mindose"
                  ,label = "Min dose"
                  ,value = "50") # ""
        , textInput(inputId = "maxdose"
                    ,label = "Max dose"
                    ,value = "300") # ""
        , textInput(inputId = "dose_units"
                    ,label = "Dose units"
                    ,value = "mg")
        , cellWidths = c("35%","35%","30%")
      ),
      splitLayout(
        numericInput(inputId = "num_doses"
                     ,label = "Number of doses"
                     ,value = 5
                     ,min = 3
                     ,max = 7
                     ,step = 1)
        , radioButtons(inputId = "range_scaling"
                       ,label = "Dose level sequence"
                       ,choices = c("geom","arith","custom")
                       ,selected = "geom"
                       ,inline = TRUE)
        , cellWidths = c("35%","65%")
      ),
      uiOutput("dose_levels"),
      ), # </fieldset>
      tags$fieldset(id="optimal-dose-heterogeneity",
      tags$legend("optimal-dose heterogeneity"),
      splitLayout(
        textInput(inputId = "median_mtd"
                  ,label = HTML("median MTD<sub>i</sub>")
                  ,value = "200") # ""
        , numericInput(inputId = "sigma_median"
                       ,label = HTML("&sigma;<sub>median</sub> (%)")
                       ,value = 50
                       ,min = 10
                       ,max = 200
                       ,step = 10)
        , numericInput(inputId = "sigma_CV"
                       ,label = HTML("&sigma;<sub>CV</sub> (%)")
                       ,value = 50
                       ,min = 10
                       ,max = 200
                       ,step = 10)
        , cellWidths = c("35%","35%","30%")
      )
      ), # </fieldset>
      splitLayout(
      actionButton("resample", label = "RESAMPLE"
                 , style = "color: #fff; background-color: #00aa00"),
      ## centered button
      div(class="flexcontainer",

          ## action button
          actionButton(inputId="startHelp", label="Tutorial", class="btn-info")
      ),
      cellWidths = c("50%","50%"),
      cellArgs = list(style = "padding: 4px", align = "center")
      ), # </splitLayout>
      tags$fieldset(id="dose-escalation-design",
      tags$legend("dose-escalation design"),
      splitLayout(
        radioButtons(inputId = "design"
                     ,label = "Method"
                     ,choices = c("3 + 3","BOIN","CRM")
                     ,selected = "3 + 3"
                     ,inline = FALSE)
        , verticalLayout(
            sliderInput(inputId = "ttr"
                       ,label = "Target toxicity rate"
                       ,min = 15
                       ,max = 33
                       ,value = 25
                       ,post = "%")
          , hidden(checkboxInput("editing_skeleton", "[hidden]", value = FALSE))
          , textOutput("J")
          )
        , verticalLayout(
            radioButtons(inputId = "cohort_size"
                        ,label = HTML("Cohort size")
                        ,choices = c("2", "3")
                        ,selected = "3"
                        ,inline = TRUE)
          , numericInput(inputId = "maxcohs_perdose"
                        ,label = HTML("Max cohs/dose")
                        ,value = 2
                        ,min = 3
                        ,max = 6
                        ,step = 1)
          , disabled(
              numericInput(inputId = "enroll_max"
                          ,label = HTML("Max enrollment")
                          ,value = NA
                          ,min = 24
                          ,max = 24
                          ,step = 1)
            )
          )
        , cellWidths = c("20%","55%","25%")
      )
      ), # </fieldset>
      tags$fieldset(id="crm-skeleton",
      tags$legend("crm skeleton"),
      uiOutput("crm_skeleton"),
      ), # </fieldset>
      tags$fieldset(id="therapeutic-index", style="padding-left: 6px; padding-right: 6px",
      tags$legend("therapeutic index"),
      sliderInput(inputId = "r0"
                  ,label = HTML("Therapeutic Index r<sub>0</sub>")
                  ,min = 1.1
                  ,max = 5.0
                  ,value = 1.5),
      ), # </fieldset>
      tags$div(class = "header", checked = NA,
               tags$p(tags$b("See:"),
                      "Norris DC. Retrospective analysis of a fatal dose-finding trial",
                      tags$a(href="https://arxiv.org/abs/2004.12755",
                             tags$i("arXiv:2004.12755 [stat.ME]")),
                      "April 2020."
               ),
               tags$p(tags$b("See also:"),
                      tags$a(href="https://dcnorris.github.io/precautionary/index.html",
                             "R package 'precautionary'")
               )
      )
    ),

    ## Show a plot of the generated distribution
    mainPanel(
      plotOutput("hyperprior"),
      htmlOutput("safety")
    )
  ),

  ## Define message handlers and overlay help-system data on now-defined UI
  tags$body(tags$script(src="help.js")),

)

server <- function(input, output, session) {

  observeEvent(input$startHelp,{
    session$sendCustomMessage(type = 'startHelp', message = list(""))
  })

  num_doses <- reactive({
    in_range <- input$num_doses %in% 3:7
    shinyFeedback::feedbackWarning("num_doses", !in_range
                                 , "Should be from 3 to 7")
    req(in_range, cancel=TRUE)
    as.integer(input$num_doses)
  })

  observeEvent(input$resample, {
    MTDi_gen()$resample()
  })

  blank_safety <- rep("--", 7)
  names(blank_safety) <- c("None", paste("Grade", 1:5), "Total")
  blank_kable <- blank_safety %>% t %>% knitr::kable(align='r') %>%
    kable_styling(position = "left", full_width = FALSE) %>%
    add_header_above(c("Expected counts per toxicity grade"=6, " "=1))

  updateManagedNumericInput <- function(session, inputId, value, ...) {
    do.call(inputId, list(value), ...)
    updateNumericInput(session, inputId = inputId, value = value)
  }

  observeEvent(input$design, {
    if (input$design == "3 + 3") {
      shinyjs::disable("ttr")
      shinyjs::disable("maxcohs_perdose")
      shinyjs::disable("cohort_size")
      shinyjs::disable("crm_skeleton")
      updateNumericInput(session, inputId = "cohort_size", value = 3)
      updateManagedNumericInput(session, "maxcohs_perdose", 2L)
      updateNumericInput(session, inputId = "enroll_max", value = NA)
    } else if (input$design == "BOIN") {
      updateManagedNumericInput(session, "maxcohs_perdose", 3L)
      updateNumericInput(session, inputId = "enroll_max", value = ENROLL_MAX)
      shinyjs::enable("ttr")
      shinyjs::enable("maxcohs_perdose")
      shinyjs::enable("cohort_size")
      shinyjs::disable("crm_skeleton")
    } else if (input$design == "CRM") {
      updateManagedNumericInput(session, "maxcohs_perdose", 3L)
      updateNumericInput(session, inputId = "enroll_max", value = ENROLL_MAX)
      shinyjs::enable("ttr")
      shinyjs::enable("maxcohs_perdose")
      shinyjs::enable("cohort_size")
      shinyjs::enable("crm_skeleton")
    } else
      stop() # all cases exhausted
  })

  dose_counter <- reactive(seq_len(num_doses())) # c(1,...,K) for K doses

  ## TODO: Perform validation of these inputs as per
  ## https://mastering-shiny.org/action-feedback.html#validate
  mindose <- reactive({
    mindose <- as.numeric(input$mindose)
    isnum <- !is.na(mindose)
    ispos <- mindose > 0
    shinyFeedback::feedbackWarning("mindose", !isnum, "Invalid dose")
    shinyFeedback::feedbackWarning("mindose", !ispos, "Dose must be > 0")
    req(isnum && ispos)
    mindose
  })
  maxdose <- reactive({
    maxdose <- as.numeric(input$maxdose)
    isnum <- !is.na(maxdose)
    gtmin <- maxdose > mindose()
    shinyFeedback::feedbackWarning("maxdose", !isnum, "Invalid dose")
    shinyFeedback::feedbackWarning("maxdose", !gtmin, "max < min!")
    req(isnum && gtmin)
    maxdose
  })

  median_mtd <- reactive({
    median_mtd <- as.numeric(input$median_mtd)
    isnum <- !is.na(median_mtd)
    ispos <- median_mtd > 0
    shinyFeedback::feedbackWarning("median_mtd", !isnum || !ispos, "Invalid dose")
    req(isnum && ispos)
    median_mtd
  })

  sigma_median <- reactive({
    sigma_median <- as.numeric(input$sigma_median)
    isnum <- !is.na(sigma_median)
    nonneg <- sigma_median >= 0
    shinyFeedback::feedbackWarning("sigma_median", !isnum || !nonneg, "Invalid")
    req(isnum && nonneg)
    sigma_median
  })

  sigma_CV <- reactive({
    sigma_CV <- as.numeric(input$sigma_CV)
    isnum <- !is.na(sigma_CV)
    nonneg <- sigma_CV >= 0
    shinyFeedback::feedbackWarning("sigma_CV", !isnum || !nonneg, "Invalid")
    req(isnum && nonneg)
    sigma_CV
  })

  ## NB: Reading mindose() & maxdose() directly avoids a *race condition*
  readDoseLevels <- reactive(
    c(mindose()
      , sapply(2:(num_doses()-1), function(k) as.numeric(input[[paste0("D",k)]]))
      , maxdose()
      )
  )

  dose_levels <- reactive(
    switch(input$range_scaling,
           geom = exp(seq(from = log(mindose())
                          , to = log(maxdose())
                          , length.out = num_doses())),
           arith = seq(from = mindose()
                       , to = maxdose()
                       , length.out = num_doses()),
           custom = readDoseLevels()
    )
  )

  output$dose_levels <- renderUI({
    do.call(splitLayout, lapply(dose_counter(), function(k, ...)
      textInput(inputId = paste0("D", k)
              , label = HTML(paste0("D<sub>", k,"</sub>"))
              , value = dose_levels()[k]
              , ...)))
  })

  cohort_size <- reactive(
    ifelse(input$design == "3 + 3", 3, as.integer(input$cohort_size))
  )

  ## This reactiveVal effectively encapsulates input$maxcohs_perdose
  ## as a 'managed input' on which I can use updateNumericInput()
  ## without spawning recalculation side-effects. Provided that I set
  ## this reactiveVal before any such update, that update (messaged,
  ## and so delayed) will amount to a NO-OP from the perspective of
  ## this reactiveVal and of its dependencies.
  maxcohs_perdose <- reactiveVal(2L) # initialized for 3+3 default
  observeEvent(input$maxcohs_perdose, {
    maxcohs_perdose(as.integer(input$maxcohs_perdose))
  })

  cohort_max <- reactive({
    if (input$design == "3 + 3")
      return(6)
    cohort_max <- maxcohs_perdose()
    toolow <- cohort_max < 3
    cohort_size <- cohort_size()
    max_enroll <- 15
    toohigh <- cohort_max * cohort_size > max_enroll
    shinyFeedback::feedbackWarning("maxcohs_perdose", toolow
                                 , "Should be &ge; 3")
    shinyFeedback::feedbackWarning("maxcohs_perdose", toohigh
                                 , paste0("Enroll/dose > ", max_enroll)
                                   )
    req(!toolow && !toohigh)
    cohort_max * cohort_size
  })

  enroll_max <- reactive(
    ifelse(input$design == "3 + 3", NA_integer_, ENROLL_MAX)
  )

  ## Like the 'managed input' maxcohs_perdose(), this crm_skeleton()
  ## serves as a dependency buffer between the input$P_ text inputs
  ## and downstream dependencies like the cpe() and hyperprior plot.
  crm_skeleton <- reactiveVal(NULL)

  observeEvent(input$editing_skeleton, {
    req(!isTruthy(input$editing_skeleton)
      , cancelOutput = TRUE) # avoid blanking the outputs
    crm_skeleton(
    sapply(1:num_doses(),
           function(k)
             as.numeric(
                 input[[paste0("P",k)]]))
    )
    MTDi_gen()$skeleton(crm_skeleton())
  })

  output$crm_skeleton <- renderUI({
    if (input$design == "CRM") {
      crm_skeleton(MTDi_gen()$skeleton())
      do.call(splitLayout, lapply(dose_counter(), function(k, ...)
        textInput(inputId = paste0("P", k)
                , label = HTML(paste0("P<sub>", k,"</sub>"))
                , value = crm_skeleton()[k]
                , ...)
        ))
    }
  })

  MTDi_gen <- reactive(
    HyperMTDi_lognormal$new(CV = 0.01*sigma_CV()
                          , median_mtd = median_mtd()
                          , median_sdlog = 0.01*sigma_median()
                          , units = input$dose_units
                          , n = 1000)$
      doses(dose_levels())$
        skeleton(isolate(crm_skeleton()))
  )

  ENROLL_MAX <- 24 # TODO: Allow some user control (easier than explaining!)

  cpe <- reactive({
    cat("\n[REACTIVE:cpe] input$design =", input$design, "\n")
    if (input$design != "CRM") # TODO: It doesn't feel entirely right to be managing
      crm_skeleton(NULL)       #       crm_skeleton() here, although it does work.
    isolate(MTDi_gen())$skeleton(crm_skeleton())
    if (input$design == "3 + 3")
      return(Cpe3_3$new(D = num_doses()))
    cohort_sizes <- rep(cohort_size(), enroll_max()/cohort_size())
    switch(input$design
         , CRM = Crm$new(skeleton = crm_skeleton()
                       , target = 0.01*input$ttr)$
             stop_func(function(x) {
               y <- stop_for_excess_toxicity_empiric(x,
                                                     tox_lim = 0.01*input$ttr + 0.1,
                                                     prob_cert = 0.72,
                                                     dose = 1)
               if(y$stop){
                 x <- y
               } else {
                 x <- dtpcrm::stop_for_consensus_reached(x, req_at_mtd = cohort_max())
                 if(!x$stop)
                   x <- dtpcrm::stop_for_sample_size(x, enroll_max())
               }
               return(x)
             })$
             no_skip_esc(TRUE)$
             no_skip_deesc(FALSE)$
             global_coherent_esc(TRUE)
         , BOIN = Boin$new(target = 0.01*input$ttr
                          ,cohort_max = cohort_max()
                          ,enroll_max = enroll_max())$
             max_dose(num_doses())
           ) -> model
    plan(multicore)
    progressr::withProgressShiny(
      model$trace_paths(root_dose = 1
                      , cohort_sizes = cohort_sizes
                      , future.scheduling = TRUE
                        )
    , message = "Running CPE..."
    )
  })

  output$J <- renderText({
    paste(cpe()$J(), "paths")
  })

  safety <- reactive({
    input$resample # take dependency
    MTDi_gen()$fractionate(cpe = cpe()
                          ,kappa = log(input$r0))
  })

  output$safety <- renderText({
    ifelse(is.null(safety()),
           blank_kable,
           safety_kable(safety()))
  })

  ## Manage the dose-level inputs
  observeEvent({
    dose_levels()
  }, {
    shinyjs::delay(0, { # https://github.com/daattali/shinyjs/issues/54#issuecomment-235347072
      shinyjs::disable("D1")
      if (input$range_scaling == "custom")
        for(k in 2:(num_doses()-1)) shinyjs::enable(paste0("D", k))
      else
        for(k in 2:(num_doses()-1)) shinyjs::disable(paste0("D", k))
      shinyjs::disable(paste0("D", num_doses()))
    })
  })

  output$hyperprior <- renderPlot({
    input$resample # take dependency
    crm_skeleton() # depend on this to toggle or move skeleton points
    log_median <- log(median_mtd())
    half_width <- 3 * 0.01*(sigma_median() + sigma_CV())
    xlim <- exp(log_median + half_width * c(-1, 1))
    MTDi_gen()$plot(col=adjustcolor("red", alpha=0.25), xlim = xlim)
  })

}

## Run the application
shinyApp(ui = ui, server = server)

