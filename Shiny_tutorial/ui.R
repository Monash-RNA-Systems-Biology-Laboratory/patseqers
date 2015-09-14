library(shiny)
library(ggplot2)
source("helper.R")


shinyUI(fluidPage(
    titlePanel("My first app"),
    
    sidebarLayout(
        sidebarPanel(

            sliderInput("num", label = h3("Number of points"), value = 1, min = 0, max = 5000)     
            ),
        
        mainPanel(
            plotOutput("example.plot")
            )
        
        )
    
    ))