library("ggplot2")
source("helper.R")


shinyServer(function(input, output) {
    
    call_table <- reactive({
        make_table(input$slider1)
    })
    
    output$plot1 <- renderPlot({
        df <- call_table()
        
        ggplot(data=df, aes(x=x, y=y)) + 
            geom_point() +
            ggtitle("Awesome Plot")
        
    })
    
    output$table1 <- renderTable({
        head(call_table(), input$table.num)
    })
    
})