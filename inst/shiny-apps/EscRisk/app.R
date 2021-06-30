##
## This was a Shiny web application providing a simple interface to the
## simulation capabilities of the 'precautionary' package.
## Now I want to obtain an 'MCSE-free' version that generates
## safety schematics such as Fig 3 of Norris 2020c.
##
## TODO:
## MARK II -- eliminate trial-wise discrete-event sim
## 1. Factor out 'escalation' in favor of CPE-based computations
##    - Use easy defaults for args not yet supported in the UI
##    - Ideally, refactor simulate_trials to use 'exact' matrix math;
##      this leaves the UI essentially unchanged!
## 2. Run longer CPEs with pop-up progress bar?
## 3. Add & activate 'consensus' enrollment to CRM/BOIN
## 4. Account properly for whatever MCSE remains
##    - Does the existing MCSE calculation carry forward
##      validly even with CPE?
##    - Does the MCSE progress bar even remain necessary,
##      or does hyperprior sampling loop run so fast that
##      I can drop it for practical purposes?
## 5. Add a Bernoulli dose scaling
##    - Find a decent *citation* for this, beyond rumor & innuendo!
## 6. Visual feedback on hyperprior sampling via top-right plot?
##    - This could be profoundly helpful in conveying the meaning
##      of what's going on under the hood! Each trace added to the
##      plot during sampling represents another 'scenario' of 'truth'
##      that is being explored and tallied.
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
                       ,max = 200)
        , numericInput(inputId = "sigma_CV"
                       ,label = HTML("&sigma;<sub>CV</sub> (%)")
                       ,value = 50
                       ,min = 10
                       ,max = 200)
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
          , textOutput("J")
          )
        , verticalLayout(
            numericInput(inputId = "cohort_size"
                        ,label = HTML("Cohort size")
                        ,value = 3
                        ,min = 2
                        ,max = 3
                        ,step = 1)
          , numericInput(inputId = "cohort_max"
                        ,label = HTML("Max / dose")
                        ,value = 6
                        ,min = 6
                        ,max = 15
                        ,step = 3)
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
                      tags$a(href="https://cran.r-project.org/web/packages/precautionary/vignettes/Intro.html",
                             "Introduction to package 'precautionary'")
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
    validate(
      need(input$num_doses %in% 3:7, "Should be an integer from 3 to 7")
    )
    as.integer(input$num_doses)
  })

  ## Like an invisible actionButton -- just bump count to recompute
  state <- reactiveValues(cpe_count = 0)

  observeEvent(input$resample, {
    cat("Resampling..")
    MTDi_gen()$resample()
    cat("..done.\n")
  })

  blank_safety <- rep("--", 7)
  names(blank_safety) <- c("None", paste("Grade", 1:5), "Total")
  blank_kable <- blank_safety %>% t %>% knitr::kable(align='r') %>%
    kable_styling(position = "left", full_width = FALSE) %>%
    add_header_above(c("Expected counts per toxicity grade"=6, " "=1))

  observeEvent(input$design, {
    if (input$design == "3 + 3") {
      shinyjs::disable("ttr")
      shinyjs::disable("cohort_max")
      shinyjs::disable("cohort_size")
      shinyjs::disable("crm_skeleton")
      updateNumericInput(session, inputId = "cohort_size", value = 3)
      updateNumericInput(session, inputId = "cohort_max", value = 6)
    } else if (input$design == "BOIN") {
      shinyjs::enable("ttr")
      shinyjs::enable("cohort_max")
      shinyjs::enable("cohort_size")
      shinyjs::disable("crm_skeleton")
    } else if (input$design == "CRM") {
      shinyjs::enable("ttr")
      shinyjs::enable("cohort_max")
      shinyjs::enable("cohort_size")
      shinyjs::enable("crm_skeleton")
      ## Also here set the skeleton from default when I first select (x) CRM
      probs <- MTDi_gen()$avg_tox_probs()
      for (k in dose_counter())
        updateTextInput(session, inputId = paste0("P",k), value = probs[k])
    } else
      stop() # all cases exhausted
  })

  observeEvent(input$cohort_size, {
    ## Actively manage the 'cohort_max' input to maintain consistency with cohort_size..
    cohort_size <- input$cohort_size
    cohort_max <- round(input$cohort_max / cohort_size) * cohort_size
    if (input$cohort_max == cohort_max) { # If update will not change value ...
      ## then the updateNumericInput below will NOT automatically trigger
      ## the input$cohort_max event below, and so we require a 'manual'
      ## triggering of design() recalculation...
      state$cpe_count <- state$cpe_count + 1
    }
    updateNumericInput(session, inputId = "cohort_max"
                     , min = c(NA, 8, 9)[cohort_size]
                     , max = c(NA, 16, 15)[cohort_size]
                     , step = cohort_size
                     , value = cohort_max)
  })

  ## Explicitly trigger design() recalculation when input$cohort_max changes,
  ## effectively 'channeling' design()'s several dependencies through this variable.
  observeEvent(input$cohort_max, {
    state$cpe_count <- state$cpe_count + 1
  })

  dose_counter <- reactive(seq_len(num_doses())) # c(1,...,K) for K doses

  ## TODO: Perform validation of these inputs as per
  ## https://mastering-shiny.org/action-feedback.html#validate
  mindose <- reactive({
    mindose <- as.numeric(input$mindose)
    isnum <- !is.na(mindose)
    ispos <- mindose > 0
    shinyFeedback::feedbackWarning("mindose", !isnum, "Invalid dose")
    shinyFeedback::feedbackWarning("mindose", !ispos, "Be serious!")
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

  ## NB: Reading mindose() & maxdose() directly avoids a *race condition*
  readDoseLevels <- reactive(
    ## TODO: A problem seems to be that this does not take a dependency
    ##       on any of the intermediate doses.
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

  crm_skeleton <- reactive(
    sapply(1:num_doses(), function(k) as.numeric(input[[paste0("P",k)]]))
  )

  output$crm_skeleton <- renderUI({
    skeleton <- isolate(crm_skeleton())
    do.call(splitLayout, lapply(dose_counter(), function(k, ...)
      if (input$design == "CRM")
        textInput(inputId = paste0("P", k)
                , label = HTML(paste0("P<sub>", k,"</sub>"))
                , value = skeleton[k]
                , ...)
      else
        disabled(
          textInput(inputId = paste0("P", k)
                  , label = HTML(paste0("P<sub>", k,"</sub>"))
                  , value = skeleton[k]
                  , ...)
        )
      ))
  })

  ## Changes to MTDi_gen don't propagate, so let me create a simple
  ## reactiveVal in the hope that solves the problem...
  ## TODO: Consider including this as component of 'state'.
  ## TODO: Find a more natural idiom.
  nsamples <- reactiveVal(1000)

  MTDi_gen <- reactive(
    HyperMTDi_lognormal$new(CV = 0.01*as.numeric(input$sigma_CV)
                          , median_mtd = as.numeric(input$median_mtd)
                          , median_sdlog = 0.01*as.numeric(input$sigma_median)
                          , units = input$dose_units
                          , n = nsamples())$doses(dose_levels())
  )

  ENROLL_MAX <- 24 # TODO: Allow some user control (easier than explaining!)

  ## Converting this expensive recalculation to 'eventReactive',
  ## in the hope this creates opportunity to explicitly trigger
  ## recomputations.
  design <- eventReactive({ state$cpe_count; input$design; dose_levels() }, {
    cohort_size <-input$cohort_size
    cat("recalculating design() with parameters:\n")
    cat("  cohort_size =", cohort_size,
        "\n  cohort_max =", input$cohort_max,
        "\n  num_doses =", num_doses(),
        "\n  dose_levels =", dose_levels(),
        "...")
    ## Okay, NOW we can proceed ...
    des <- switch(input$design,
           `3 + 3` = list(b = precautionary:::b[[num_doses()]]
                         ,U = precautionary:::U[[num_doses()]]
                          )
           ## TODO: Provide UI inputs for the CRM skeleton, and display on plot
         , CRM = Crm$new(skeleton = MTDi_gen()$avg_tox_probs()
                       , target = 0.01*input$ttr)$
             stop_func(function(x) {
               y <- stop_for_excess_toxicity_empiric(x,
                                                     tox_lim = 0.01*input$ttr + 0.1,
                                                     prob_cert = 0.72,
                                                     dose = 1)
               if(y$stop){
                 x <- y
               } else {
                 x <- dtpcrm::stop_for_consensus_reached(x, req_at_mtd = input$cohort_max)
                 if(!x$stop)
                   x <- dtpcrm::stop_for_sample_size(x, ENROLL_MAX)
               }
               return(x)
             })$
             no_skip_esc(TRUE)$
             no_skip_deesc(FALSE)$
             global_coherent_esc(TRUE)$
             trace_paths(root_dose=1
                       , cohort_sizes=rep(cohort_size, ENROLL_MAX/cohort_size)
                       , impl = 'rusti')

         , BOIN = Boin$new(target = 0.01*input$ttr
                          ,cohort_max = input$cohort_max
                          ,enroll_max = ENROLL_MAX)$
             max_dose(num_doses())$
             trace_paths(root_dose=1
                       , cohort_sizes=rep(cohort_size, ENROLL_MAX/cohort_size)
                         )

    ) # </switch>
    cat("\n")
    return(des)
  })

  output$J <- renderText({
    cpe <- design()
    if (is(cpe,'Cpe')) {
      J <- cpe$J()
    } else {
      J <- length(cpe$b)
    }
    paste(J, "paths")
  })

  safety <- reactive({
    input$resample # take dependency
    cpe <- design()
    if (is(cpe, 'Cpe'))
      cpe <- cpe$bU()
    MTDi_gen()$fractionate(b = cpe$b
                          ,U = cpe$U
                          ,kappa = log(input$r0))
  })

  output$safety <- renderText({
    input$resample # take dependency
    ifelse(is.null(safety()),
           blank_kable,
           safety_kable(safety()))
  })

  ## Any one of these many UI events will invalidate the safety table:
  observeEvent({
    dose_levels()
    MTDi_gen() # TODO: Does this even generate events, presently?
    nsamples()
    design()
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

  ## TODO: A nicety worth considering is to keep the horizontal axis fixed
  ##       between resamplings. This helps the user to see more clearly
  ##       how multiple resamplings compare. Presently, the axis recomputes
  ##       with each resampling, jolting the plot misleadingly.
  output$hyperprior <- renderPlot({
    input$resample # take dependency
    MTDi_gen()$plot(col=adjustcolor("red", alpha=0.25))
  })

}

## Run the application
shinyApp(ui = ui, server = server)

