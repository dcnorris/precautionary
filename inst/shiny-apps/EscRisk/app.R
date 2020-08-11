#
# This is a Shiny web application providing a simple interface to the
# simulation capabilities of the 'precautionary' package.
#
# TODO:
#/1. Inputs for min, max, # doses, geometric vs arithmetic sequence
#/2. Inputs for median(MTDi), sigma_med and sigma_CV.
# 3. Simulate 3+3 trial for M=100 reps, show table
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

library(shiny)
library(precautionary)

ui <- fluidPage(
  
  includeCSS("www/tweaks.css"),
  
  # Application title
  titlePanel("Predict Risks of High-Grade Toxicities in Dose-Escalation"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      splitLayout(
        textInput(inputId = "mindose"
                  ,label = "Min. dose"
                  ,value = "")
        , textInput(inputId = "maxdose"
                    ,label = "Max. dose"
                    ,value = "")
        , textInput(inputId = "dose_units"
                    ,label = "Dose unit"
                    ,value = "mg")
        , cellWidths = c("35%","35%","30%")
      ),
      splitLayout(
        numericInput(inputId = "num_doses"
                     ,label = "Number of doses"
                     ,value = 4
                     ,min = 2
                     ,max = 8)
        , radioButtons(inputId = "range_scaling"
                       ,label = "Dose level sequence"
                       ,choices = c("geometric","arithmetic")
                       ,selected = "geometric"
                       ,inline = TRUE)
      ),
      hr(),
      splitLayout(
        textInput(inputId = "median_MTD"
                  ,label = HTML("median MTD<sub>i</sub>")
                  ,value = "")
        , numericInput(inputId = "sigma_med"
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
        , sliderInput(inputId = "target_tox_rate"
                      ,label = "Target toxicity rate for CRM or BOIN (%)"
                      ,min = 15
                      ,max = 33
                      ,value = 25)
        , cellWidths = c("40%","60%")
      ),
      sliderInput("CV"
                    ,"Select coefficient of variation (CV) for MTDi in population"
                    ,min = 0.2
                    ,max = 1.4
                    ,value = 0.75),
      sliderInput("theMTD"
                  ,"Adjust 'the' MTD to maximize population-level efficacy"
                  ,min = 50
                  ,max = 500
                  ,value = 200),
      checkboxInput("rainbow"
                    ,"Show MTDi rainbow"
                    ,value = FALSE),
      hr(),
      tags$div(class = "header", checked = NA,
               tags$p(tags$b("See:"), "Norris DC. Costing 'the' MTD.",
                      tags$i("bioRxiv"), "August 2017:150821.",
                      tags$a(href="http://www.biorxiv.org/content/early/2017/08/22/150821",
                             "doi:10.1101/150821.")
               )
      )
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot")
    )
  )
)

server <- function(input, output) {
  
  mean.MTDi <- 200 # i.e., alpha/beta = shape*scale
  shape <- reactive((input$CV)^-2)
  scale <- reactive(mean.MTDi/input$theMTD/shape())
  # Precompute the density curve each time the input$CV changes.
  # This precomputed density needs to span the full width of theta plotted,
  # at the lowest MTDthe on the slider.
  density <- reactive(dinvgamma((0:2000)/1000
                                , shape = (input$CV)^-2
                                , scale = mean.MTDi/50/(input$CV)^-2
  ))
  # That precomputed 'density' should now be indexable by the following function:
  densityFun <- function(theta){
    (50/input$theMTD)*density()[round(1000 * theta*50/input$theMTD) + 1]
  }
  Pr <- reactive({
    alpha <- shape()
    beta <- 1/scale()
    Q <- Rgamma(alpha - 0.5, beta, lower=FALSE) # Regularized gamma function
    0.5 * sqrt(beta)*gamma(alpha - 0.5)/gamma(alpha) * Q
  })
  
  # Build non-reactive data.frame with all series for plotting
  x <- (0:200)/100
  x. <- c(x[x<=1], x[x>=1])
  X <- x[x<=1]
  Z <- rep(0.0, 1+length(x[x>1]))
  Emax <- rbind(data.frame(x=X, Emax=0.5*sqrt(X), series="approx"),
                data.frame(x=x[x>1], Emax=(x/(x+1))[x>1], series="Emax"))
  Emax$series <- as.factor(Emax$series)
  
  output$distPlot <- renderPlot({
    right <- xYplot(Emax ~ x, group=series, data=Emax, type="l", lty=c(1,2)
                    , xlim=c(0,2), ylim=c(0,0.7)
                    , ylab = "Probability of achieving remission"
                    , label.curves = FALSE
                    , panel = function(x,y,...) {
                      panel.points(x=1, y=0.5, pch=10, cex=1.2)
                      panel.arrows(x0=1, y0=0.5, x1=2, y1=0.5, lty=1, type="closed", length=0.1)
                      panel.text(x=2, y=0.52, adj=c(1,NA),
                                 "Remission rate for individualized 'MTDi' dosing  ")
                      panel.arrows(x0=0, y0=Pr(), x1=2, y1=Pr(), lty=3, type="closed", length=0.1)
                      panel.text(x=2, y=Pr()+0.02, adj=c(1,NA),
                                 paste("Population-level remission rate under 1-size-fits-all dosing at 'the' MTD =",
                                       input$theMTD, " "))
                      panel.arrows(x0=1.25, y0=0.48, x1=1.25, y1=Pr()+0.05, lty=1, lwd=2, type="closed",
                                   length=0.15)
                      panel.text(x=1.27, y=(Pr()+0.5)/2, adj=c(0,0),
                                 paste(round(200*(0.5-Pr())), "% loss of population-level efficacy", sep=''))
                      panel.xYplot(x,y,...)
                    })
    # 2. Overplot Inv-Gamma dashed with chopped version
    #shape <- (input$CV)^-2
    #scale <- mean.MTDi/input$theMTD/shape()
    sub.half <- Rgamma(a=shape(), x=(0.5*scale())^-1, lower=FALSE)
    excluded <- Rgamma(a=shape(), x=(1*scale())^-1, lower=TRUE)
    dist <- rbind(data.frame(x=x, f=densityFun(x)))
    dist$rainbow <- rev(rainbow(n=nrow(dist), start=0, end=0.7, alpha=0.35))
    dissed <- rbind(data.frame(x=x., f=c(densityFun(X), Z)))
    
    left <- xYplot(f ~ x, data=dist, type='l', lty=2
                   , ylim=c(0,2.4) # TODO: Recalculate to avoid clipping?
                   , ylab="Density"
                   , panel=function(x,y,...,cutoff){
                     m1 <- min(which(x >= 0.0))
                     m2 <- max(which(x <= cutoff))
                     tmp <- data.frame(x1 = x[c(m1,m1:m2,m2)], y1 = c(0,y[m1:m2],0))
                     panel.polygon(tmp$x1, tmp$y1, col="lightgray", border="lightgray")
                     tailTextCol <- ifelse(input$rainbow
                                           , 'black'
                                           , trellis.par.get('superpose.line')$col)
                     panel.text(x=0.35, y=0.2
                                , col = tailTextCol
                                , paste(round(100*sub.half),
                                        "% treated\nat under half\ntheir MTDi", sep=''))
                     panel.text(x=1.23, y=0.2
                                , col = tailTextCol
                                , paste(round(100*excluded),
                                        "% untreated\nbecause 'the' MTD\nexceeds their MTDi", sep=''))
                     # TODO: Draw the length(x)-1 polygons in various colors
                     if(input$rainbow){
                       for(i in 1:length(x[-1])){
                         # Compute cumulative frequency as index into rainbow
                         F <- pinvgamma(i/length(x), shape=shape(), scale=scale())
                         with(dist,
                              panel.polygon(x = rep(x[i:(i+1)], each=2)
                                            ,y = c(0, f[i:(i+1)], 0)
                                            ,col = rainbow[floor(F^0.4*length(x))+1]
                                            ,lty = 0
                              )
                         )
                       }
                     }
                     panel.xYplot(x,y,...)
                   }
                   , cutoff=0.5 # shade all treated at under half their MTDi's
    )
    left <- left + xYplot(f ~ x, data=dissed, type='l', lty=1)
    update(doubleYScale(left, right, add.ylab2=TRUE)
           , par.settings = simpleTheme(col=brewer.pal(4,"PRGn")[c(4,1)])
           , xlab=expression(theta[i] %==% MTD[the] / MTD[i])
           , xlim=c(0,2.0)
    )
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

