# This code shows how the UI will be laid out
shinyUI(fluidPage(
  titlePanel("Making some points"),
  # Inputs
  sidebarLayout(
    sidebarPanel(
      
      sliderInput("slider1", label = h3("Show N points"), 
                  min = 0, max = 500, value = 1),
      
      numericInput("num", label = h3("Show N rows"),
                   value = 1, max = 500)
    ),
    # Outputs
    mainPanel(
      plotOutput("plot1"),
      tableOutput("table1")
    )
  )
  
))