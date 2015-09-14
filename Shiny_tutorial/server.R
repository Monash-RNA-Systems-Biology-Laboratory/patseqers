library(shiny)
library(ggplot2)
source("helper.R")


shinyServer(function(input, output) {
    
    selection.x <- reactive({
        return(input$Xaxis)
    })
    
    selection.y <- reactive({
        return(input$Yaxis)
    })
#     
    output$example.plot <- renderPlot({
        selX <- selection.x()
        selY <- selection.y()

        ggplot(df, aes_string(x = selX, y = selY)) + geom_point()
    })
    
    output$plot.info <- renderText({
        paste0("x=", input$plot_click$x, "y=", input$plot_click$y)
    })
    
})