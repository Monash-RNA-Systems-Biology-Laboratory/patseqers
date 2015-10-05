# Loads everything in helper.R into memory 
source("helper.R")
# The basic setup of the shiny server function
shinyServer(function(input, output) {

  df_values <- reactive({    
    df <- make_df(input$slider1)      
  })

  output$plot1 <- renderPlot({ 
    df <- df_values()
    return(plot(df$x,df$y, main = "Do you see my Points?"))          
  })

  output$table1 <- renderTable({
    df <- df_values()
    tab <- head(df, input$num)
    return(tab)
  })
  
})





