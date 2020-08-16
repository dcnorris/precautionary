#
# This is a Shiny web application providing a simple interface to the
# simulation capabilities of the 'precautionary' package.
#
# TODO:
#x1. Logarithmic scaling of TI slider
#  -> Remarkably, no such slider is available off-the-shelf!
# 2. Craft a complete intro/tutorial https://shiny.rstudio.com/articles/js-introjs.html
# (a) Implement the intro using package::introJS
# (b) Factor out the dependency on this (non-CRAN!) package

library(shiny)
library(precautionary)
library(kableExtra)

ui <- fluidPage(
  shinyjs::useShinyjs(),
  shinyFeedback::useShinyFeedback(),

  includeCSS("www/introjs.css"), # TODO: includeCSS("www/introjs.min.css"),
  includeScript("www/intro.js"), # TODO: includeScript("www/intro.min.js"),
  
  # Message handlers and help-system data
  includeScript("www/help.js"),

  includeCSS("www/tweaks.css"),
  
  # Application title
  titlePanel("Predict Risks of High-Grade Toxicities in Dose-Escalation Trials"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      splitLayout(
        textInput(inputId = "mindose"
                  ,label = "Min. dose"
                  ,value = "50") # ""
        , textInput(inputId = "maxdose"
                    ,label = "Max. dose"
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
                     ,min = 2
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
      hr(),
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
      ),
      hr(),
      splitLayout(
        radioButtons(inputId = "design"
                     ,label = "Dose escalation"
                     ,choices = c("3 + 3","CRM","BOIN")
                     ,selected = "3 + 3"
                     ,inline = FALSE)
        , sliderInput(inputId = "ttr"
                      ,label = "Target toxicity rate"
                      ,min = 15
                      ,max = 33
                      ,value = 25
                      ,post = "%")
        , numericInput(inputId = "enroll_max"
                       ,label = HTML("Max. enroll")
                       ,value = 24
                       ,min = 18
                       ,max = 36
                       ,step = 3)
        , cellWidths = c("30%","45%","25%")
      ),
      uiOutput('RunStopButton'),
      # centered button
      div(class="flexcontainer", 
          
          # action button
          actionButton(inputId="startHelp", label="start", class="btn-success")
      ),
      hr(),
      sliderInput(inputId = "r0"
                  ,label = HTML("Therapeutic Index r<sub>0</sub>")
                  ,min = 1.1
                  ,max = 5.0
                  ,value = 1.5),
      hr(),
      tags$div(class = "header", checked = NA,
               tags$p(tags$b("See:"),
                      "Norris DC. Retrospective analysis of a fatal dose-finding trial",
                      tags$a(href="http://arxiv.org/abs/2004.12755",
                             tags$i("arXiv:2004.12755 [stat, q-bio]")),
                      "April 2020."
               ),
               tags$p(tags$b("See also:"),
                      tags$a(href="https://cran.r-project.org/web/packages/precautionary/vignettes/Intro.html",
                             "Introduction to package 'precautionary'")
               )
      )
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("hyperprior"),
      plotOutput("simprogress", height = "150px"),
      htmlOutput("safety")
    )
  )
)

server <- function(input, output, session) {
  
  observeEvent(input$startHelp,{
    session$sendCustomMessage(type = 'startHelp', message = list(""))
  })
  
  # Let me try implementing a self-toggling Run/Stop button in 'pure Shiny',
  # without exploiting Javascript. This would seem to require maintaining the
  # desired state as a reactiveVal, and re-rendering the button accordingly.
  state <- reactiveValues(sim = 'ready') # 'ready' | 'running'
  
  output$RunStopButton <- renderUI({
    if (state$sim == 'ready')
      actionButton("run_stop", label = "RUN simulation"
                   , style = "color: #fff; background-color: #00aa00")
    else
      actionButton("run_stop", label = "STOP simulation"
                   , style = "color: #fff; background-color: #cd0000")
  })
  
  observeEvent(input$run_stop, {
    if (isolate(state$sim) == 'ready')
      state$sim <- 'running'
    else
      state$sim <- 'ready'
  })
  
  blank_safety <- rep("--", 7)
  names(blank_safety) <- c("None", paste("Grade", 1:5), "Total")
  blank_kable <- blank_safety %>% t %>% knitr::kable(align='r') %>%
    kable_styling(position = "left", full_width = FALSE) %>%
    add_header_above(c("Expected counts per toxicity grade"=6, " "=1))
  
  options(ordinalizer = function(dose, r0) { # NB: no r0 default
    c(`Grade 1` = dose/r0^2
      , `Grade 2` = dose/r0
      , `Grade 3` = dose
      , `Grade 4` = dose*r0
      , `Grade 5` = dose*r0^2)
  })
  
  observeEvent(input$design, {
    if (input$design == "3 + 3") {
      shinyjs::disable("ttr")
      shinyjs::disable("enroll_max")
    } else {
      shinyjs::enable("ttr")
      shinyjs::enable("enroll_max")
    }
  })
  
  dose_counter <- reactive(seq_len(input$num_doses)) # c(1,...,K) for K doses
  
  # TODO: Perform validation of these inputs as per
  # https://mastering-shiny.org/action-feedback.html#validate
  mindose <- reactive({
    mindose <- as.numeric(input$mindose)
    isnum <- !is.na(mindose)
    ispos <- mindose > 0
    shinyFeedback::feedbackWarning("mindose", !isnum, "Invalid dose")
    shinyFeedback::feedbackWarning("mindose", !ispos, "Seriously?")
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
  
  # NB: Reading mindose() & maxdose() directly avoids a *race condition*
  readDoseLevels <- reactive(
    c(mindose()
      , sapply(2:(input$num_doses-1), function(k) as.numeric(input[[paste0("D",k)]]))
      , maxdose()
      )
  )
  
  dose_levels <- reactive(
    switch(input$range_scaling,
           geom = exp(seq(from = log(mindose())
                          , to = log(maxdose())
                          , length.out = input$num_doses)),
           arith = seq(from = mindose()
                       , to = maxdose()
                       , length.out = input$num_doses),
           custom = readDoseLevels()
    )
  )
  
  observe(
    options(dose_levels = dose_levels())
  )
  
  output$dose_levels <- renderUI({
    do.call(splitLayout, lapply(dose_counter(), function(k, ...)
      textInput(inputId = paste0("D", k)
                , label = HTML(paste0("D<sub>", k,"</sub>"))
                , value = dose_levels()[k]
                , ...)))
  })
  
  # TODO: Should the inputs be validated here?
  mtdi_gen <- reactive(
    hyper_mtdi_lognormal(CV = 0.01*as.numeric(input$sigma_CV)
                         , median_mtd = as.numeric(input$median_mtd)
                         , median_sdlog = 0.01*as.numeric(input$sigma_median)
                         , units = input$dose_units)
  )
  
  design <- reactive(
    switch(input$design,
           `3 + 3` = get_three_plus_three(num_doses = input$num_doses),
           CRM = get_dfcrm(
             skeleton = rowMeans( # TODO: Abstract this to 'precautionary' package
               sapply(precautionary:::draw_samples(mtdi_gen(), n=100)
                      , function(mtdi) mtdi@dist$cdf(dose_levels()))),
             target = 0.01*input$ttr) %>%
             stop_at_n(n = input$enroll_max),
           BOIN = get_boin(
             num_doses = 5,
             target = 0.01*input$ttr) %>%
             stop_at_n(n = 24)
    )
  )
  
  # Unlike all other expressions in this app, 'hsims' cannot be recalculated feasibly,
  # and requires incremental *updates*. This necessitates its treatment as a reactiveVal.
  hsims <- reactiveVal(NULL)
  
  safety <- reactive({
    if (is.null(hsims()))
      return(NULL)
    summary(hsims(), r0 = input$r0)$safety
  })
  
  worst_mcse <- reactive({
    if (is.null(safety()))
      return(NA)
    max(safety()['MCSE', 2:6])
  })
  
  observe({
    if (state$sim == 'running') {
      hsims(
        if (!is.null(isolate(hsims()))) {
          isolate(hsims()) %>% extend(num_sims = 10) # <11 ==> no txtProgressBar
        } else { # iteration base case
          isolate(design()) %>%
            simulate_trials(
              num_sims = 10,
              true_prob_tox = isolate(mtdi_gen())
            )
        }
      )
      invalidateLater(500, session) # repeat
    }
  })
  
  output$safety <- renderText({
    ifelse(is.null(safety()),
           blank_kable,
           safety_kable(safety()))
  })
  
  observeEvent(worst_mcse(), {
    if (!is.na(worst_mcse()) && worst_mcse() < 0.05) state$sim <- 'ready' # HALT sim
  })
  
  plotProgress <- function(worst_mcse) {
    reps_so_far <- length(hsims()$fits)
    reps_needed <- ceiling(length(hsims()$fits) * ( worst_mcse / 0.05 )^2)
    fraction_done <- reps_so_far / reps_needed
    barplot(fraction_done, width = 0.9
            , ylim = c(0,1), xlim = c(0,1)
            , horiz = TRUE, asp = 0.04
            , xlab = "Largest MCSE for expected toxicity counts, Grades 1-5"
            , main = paste0(reps_so_far, " trials")
            , axes = FALSE # will be drawn separately
    )
    tics <- c(Inf, seq(0.1, 0.05, -0.01))
    axis(side = 1, at = (0.05/tics)^2
         , labels = c(expression("" %+-% infinity), paste0("±", substring(tics[-1],2))))
  }
  
  output$simprogress <- renderPlot(plotProgress(worst_mcse()))

  # Any one of these many UI events will invalidate the safety table:
  observeEvent({
    dose_levels()
    mtdi_gen()
    design()
  }, {
    hsims(NULL)
    shinyjs::delay(0, { # https://github.com/daattali/shinyjs/issues/54#issuecomment-235347072
      shinyjs::disable("D1")
      if (input$range_scaling == "custom")
        for(k in 2:(input$num_doses-1)) shinyjs::enable(paste0("D", k))
      else
        for(k in 2:(input$num_doses-1)) shinyjs::disable(paste0("D", k))
      shinyjs::disable(paste0("D", input$num_doses))
    })
  })
  
  output$hyperprior <- renderPlot({
    set.seed(2020) # avoid distracting dance of the samples
    options(dose_levels = dose_levels())
    plot(mtdi_gen(), n=100, col=adjustcolor("red", alpha=0.25))
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)

