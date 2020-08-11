#
# This is a Shiny web application providing a simple interface to the
# simulation capabilities of the 'precautionary' package.
#
# TODO:
#/1. Inputs for min, max, # doses, geometric vs arithmetic sequence
#/2. Inputs for median(MTDi), sigma_med and sigma_CV.
#/3. Simulate 3+3 trial for M=100 reps, show table
#  -> https://shiny.rstudio.com/articles/render-table.html
# 4. Feedback on actual doses triggered by selection of geom/arith
#  -> Consider for-now disabled textInputs as means for display
# 5. Option to specify CRM or BOIN with target toxicity rate
# 6. Recursively update table until standard errors < 0.05
# 7. Show a progress bar marked in standard errors
#  -> https://shiny.rstudio.com/articles/progress.html
# 8. Use tabSetPanel to show details like hyperprior draws? (TMI?)
#/9. Improve spacing via CSS
# 10. Pop-up (or roll-over?) explanations -- via (?) or (i) symbol
# 11. Inactivate TTL when 3+3 design selected
# 12. Option to specify n doses verbatim
#  -> Could this be via enabled editing of feedback area?
# 13. Foolproof inputs constraints & checks!
# 14. Implement *reactivity* to improve clarity & performance

library(shiny)
library(precautionary)

ui <- fluidPage(
  shinyjs::useShinyjs(),
  
  includeCSS("www/tweaks.css"),
  
  # Application title
  titlePanel("Predict Risks of High-Grade Toxicities in Dose-Escalation"),
  
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
                       ,choices = c("geometric","arithmetic")
                       ,selected = "geometric"
                       ,inline = TRUE)
        , cellWidths = c("40%","60%")
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
                       ,label = "Dose-escalation design"
                       ,choices = c("3 + 3","CRM","BOIN")
                       ,selected = "3 + 3"
                       ,inline = TRUE)
        , sliderInput(inputId = "ttr"
                      ,label = "Target toxicity rate (%)"
                      ,min = 15
                      ,max = 33
                      ,value = 25)
        , cellWidths = c("50%","50%")
      ),
      splitLayout(
        actionButton(inputId = "simulate"
                     ,label = "Simulate") # TODO: Change to "Continue sim" after [Halt]
        , actionButton(inputId = "halt"
                     ,label = "Halt")
      ),
      hr(),
      sliderInput(inputId = "r0"
                  ,label = HTML("r<sub>0</sub> parameter of ordinalizer")
                  ,min = 1.1
                  ,max = 5.0
                  ,value = 1.5),
      hr(),
      tags$div(class = "header", checked = NA,
               tags$p(tags$b("See:"),
                      "Norris DC. Retrospective analysis of a fatal dose-finding trial",
                      tags$i("arXiv:2004.12755 [stat, q-bio]."),
                      "April 2020.",
                      tags$a(href="http://arxiv.org/abs/2004.12755",
                             "http://arxiv.org/abs/2004.12755")
               )
      )
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      htmlOutput("safety")
      #tableOutput("tbl")
    )
  )
)

server <- function(input, output) {
  
  observeEvent(input$design, {
    if(input$design == "3 + 3")
      shinyjs::disable("ttr")
    else
      shinyjs::enable("ttr")
  })
  
  # TODO: Try building HTML with subscripts for dose numbering: D_1, ..., D_n
  dose_display <- reactive(paste0("D", seq_len(input$num_doses)))
  
  mindose <- reactive(as.numeric(input$mindose)) # TODO: do validation within these
  maxdose <- reactive(as.numeric(input$maxdose)) #       reactive computations?
  
  dose_levels <- reactive(
    switch(input$range_scaling,
           geometric = exp(seq(from = log(mindose())
                               , to = log(maxdose())
                               , length.out = input$num_doses)),
           arithmetic = seq(from = mindose()
                            , to = maxdose()
                            , length.out = input$num_doses))
  )
  
  output$dose_levels <- renderUI({
    do.call(splitLayout, lapply(seq_len(input$num_doses), function(k, ...)
      textInput(inputId = paste0("D", k)
                , label = HTML(paste0("D<sub>", k,"</sub>"))
                , value = dose_levels()[k]
                , ...)))
  })
  
  # observeEvent(output$dose_levels, {
  #   shinyjs::disable("D1")
  # })
  
  observeEvent(input$simulate, {
    # The 'main event'!
    cat(file = stderr(), "input$mindose:", input$mindose, "\n")
    shiny::validate(
      need(0 < (mindose <- as.numeric(input$mindose))
           , "Please provide a positive minimum dose")
      ,need(mindose < (maxdose <- as.numeric(input$maxdose))
            , "Please provide a maximum dose greater than minimum")
      ,need(0 < (median_mtd <- as.numeric(input$median_mtd))
            , "Please input a valid median MTD")
      ,need(0 < (sigma_median <- 0.01*as.numeric(input$sigma_median))
            , HTML("Please input a valid &sigma;<sub>median</sub>"))
      ,need(0 < (sigma_CV <- 0.01*as.numeric(input$sigma_CV))
            , HTML("Please input a valid &sigma;<sub>CV</sub>"))
    )
    options(dose_levels = dose_levels())
    # TODO: Would separately reactive dose_levels & mtdi_gen lend clarity?
    mtdi_gen <- hyper_mtdi_lognormal(CV = sigma_CV
                                     , median_mtd = median_mtd
                                     , median_sdlog = sigma_median
                                     , units = input$dose_units)
    hsims <- get_three_plus_three(num_doses = input$num_doses) %>%
      simulate_trials(
        num_sims = 100 # ... to start
      , true_prob_tox = mtdi_gen)
    summary(hsims, ordinalizer = function(dose, r0 = input$r0) {
      c(Gr1=dose/r0^2, Gr2=dose/r0, Gr3=dose, Gr4=dose*r0, Gr5=dose*r0^2)
    })$safety -> safety
    # hsims <- hsims %>% extend(num_sims = num_sims)
    # summary(hsims, ordinalizer = function(dose, r0 = sqrt(2))
    #   c(Gr1=dose/r0^2, Gr2=dose/r0, Gr3=dose, Gr4=dose*r0, Gr5=dose*r0^2)
    # )$safety
    print(safety)
    #print(safety_kable(safety))
    #output$tbl <- renderTable(safety)
    output$safety <- renderText(safety_kable(safety))
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)

