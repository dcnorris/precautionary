# Running Shiny app 'EscRisk' ...
old <- options(device.ask.default = FALSE)
appDir <- system.file("shiny-apps", "EscRisk", package = "precautionary")
shiny::runApp(appDir, display.mode = "normal")
options(old)
