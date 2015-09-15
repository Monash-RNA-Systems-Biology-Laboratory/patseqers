shinyUI(fluidPage(
    titlePanel("My first app"),
    
    sidebarLayout(
        sidebarPanel(
            sliderInput("slider1", label = h3("Number of points"), min = 0, 
                        max = 500, value = 1),
            numericInput("table.num", label = h3("How many rows to view for the table"), value = 1, max = 500)
        ),
        
        mainPanel(
            plotOutput("plot1"),
            tableOutput("table1")
        )
    )
    
))