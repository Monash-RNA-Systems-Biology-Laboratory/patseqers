library("ggplot2")
source("helper.R")


shinyServer(function(input, output) {
    
    output$plot1 <- renderPlot({
        df <- make_table(input$slider1)
        ggplot(data= df, aes(x=x, y=y)) + geom_point() + ggtitle("Shiny app tutorial completed!")
        
    })
    
    output$table1 <- renderTable({
        df <- make_table(input$slider1)
        head(df, input$num)
    })
    
})