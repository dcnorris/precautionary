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
                          ,label = HTML("Max enroll")
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
      shinyjs::disable("maxcohs_perdose")
      shinyjs::disable("cohort_size")
      shinyjs::disable("crm_skeleton")
      updateNumericInput(session, inputId = "cohort_size", value = 3)
      updateNumericInput(session, inputId = "maxcohs_perdose", value = 2)
      updateNumericInput(session, inputId = "enroll_max", value = NA)
    } else if (input$design == "BOIN") {
      shinyjs::enable("ttr")
      shinyjs::enable("maxcohs_perdose")
      updateNumericInput(session, inputId = "maxcohs_perdose", value = 3)
      updateNumericInput(session, inputId = "enroll_max", value = ENROLL_MAX)
      shinyjs::enable("cohort_size")
      shinyjs::disable("crm_skeleton")
    } else if (input$design == "CRM") {
      shinyjs::enable("ttr")
      shinyjs::enable("maxcohs_perdose")
      updateNumericInput(session, inputId = "maxcohs_perdose", value = 3)
      updateNumericInput(session, inputId = "enroll_max", value = ENROLL_MAX)
      shinyjs::enable("cohort_size")
      shinyjs::enable("crm_skeleton")
      ## Also here set the skeleton from default when I first select (x) CRM
      probs <- MTDi_gen()$avg_tox_probs()
      for (k in dose_counter())
        updateTextInput(session, inputId = paste0("P",k), value = probs[k])
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

  cohort_max <- reactive({
    cohort_max <- as.integer(input$maxcohs_perdose)
    toolow <- cohort_max < 3
    cohort_size <- as.integer(input$cohort_size)
    max_enroll <- 15
    toohigh <- cohort_max * cohort_size > max_enroll
    shinyFeedback::feedbackWarning("maxcohs_perdose", toolow
                                 , "Should be &ge; 3")
    shinyFeedback::feedbackWarning("maxcohs_perdose", toohigh
                                 , paste0("Enroll/dose > ", max_enroll)
                                   )
    req(!toolow && !toohigh)
    cohort_max * cohort_size # TODO: Rationalize naming of input!
  })

  enroll_max <- reactive(input$enroll_max)

  crm_skeleton <- reactive({
    ## TODO: Ideally, the invalidated downstream outputs would appear
    ##       grayed to clearly announce their invalidated status. But
    ##       this is close enough for now!
    req(!isTruthy(input$editing_skeleton)
      , cancelOutput = TRUE) # avoid blanking the outputs
    sapply(1:num_doses(),
           function(k)
             as.numeric(
               isolate( # avoid over-eager recalc from a direct dependency
                 input[[paste0("P",k)]])))
  })

  output$crm_skeleton <- renderUI({
    if (input$design == "CRM") {
      auto_skeleton <- MTDi_gen()$avg_tox_probs()
      do.call(splitLayout, lapply(dose_counter(), function(k, ...)
        textInput(inputId = paste0("P", k)
                , label = HTML(paste0("P<sub>", k,"</sub>"))
                , value = auto_skeleton[k]
                , ...)
        ))
    }
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

  cpe <- reactive({
    state$cpe_count # take dependency
    cohort_size <- as.integer(input$cohort_size)
    cat(input$design, "with parameters:\n")
    cat("  cohort_size =", cohort_size,
        "\n  cohort_max =", cohort_max(),
        "\n  num_doses =", num_doses(),
        "\n  dose_levels =", dose_levels())
    if (input$design == "CRM") {
      cat("\n  skeleton =", unlist(crm_skeleton()))
    }
    cat(" ...")
    ## Okay, NOW we can proceed ...
    des <- switch(input$design,
           `3 + 3` = Cpe3_3$new(D = num_doses())
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
                 x <- dtpcrm::stop_for_consensus_reached(x, req_at_mtd = cohort_max())
                 if(!x$stop)
                   x <- dtpcrm::stop_for_sample_size(x, enroll_max())
               }
               return(x)
             })$
             no_skip_esc(TRUE)$
             no_skip_deesc(FALSE)$
             global_coherent_esc(TRUE)$
             trace_paths(root_dose=1
                       , cohort_sizes=rep(cohort_size, enroll_max()/cohort_size)
                       , impl = 'rusti')

         , BOIN = Boin$new(target = 0.01*input$ttr
                          ,cohort_max = cohort_max()
                          ,enroll_max = enroll_max())$
             max_dose(num_doses())$
             trace_paths(root_dose=1
                       , cohort_sizes=rep(cohort_size, enroll_max()/cohort_size)
                         )

    ) # </switch>
    cat("\n")
    return(des)
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

  ## Show we receive this JavaScript event
  observeEvent(input$editing_skeleton, {
    cat("skeleton being edited:", input$editing_skeleton, "\n")
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

