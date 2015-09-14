library(shiny)
library(ggplot2)
source("helper.R")


shinyUI(fluidPage(
    titlePanel("My first app"),
    
    sidebarLayout(
        sidebarPanel(

            numericInput("num", label = h3("Numeric input"), value = 1)         
            ),
        
        mainPanel(
            plotOutput("example.plot")
            )
        
        )
    
    ))