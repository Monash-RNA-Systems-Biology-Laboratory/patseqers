library(shiny)
library(ggplot2)
source("helper.R")


shinyServer(function(input, output) {
  
   
    output$example.plot <- renderPlot({
        make_plot(input$num)
    })
    
    output$plot.info <- renderText({
        paste0("x=", input$plot_click$x, "y=", input$plot_click$y)
    })
    
})