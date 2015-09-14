library(shiny)
library(ggplot2)
source("helper.R")


shinyServer(function(input, output) {
  
   
    output$example.plot <- renderPlot({
        make_plot(input$num)
    })
    
   
    
})