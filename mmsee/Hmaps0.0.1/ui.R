#ui.R
#Authors - Michael See


library(shiny)
library(DT)

shinyUI(fluidPage(
    
    title = "Heatmap",
    h1("Heatmap"),
    #uiOutput("plotui"),
    plotOutput("plot1", width="1400", height="850"),
    radioButtons("tab", label="Use Table", choices=list("False"=0, "True"=1), selected = 0, inline = T),
    
    
    hr(),
    
    conditionalPanel(condition="input.tab == 0",
                     fluidRow(
                         column(3,
                                h4("Pre-plot Information"),
                                uiOutput("selDataSet"),
                                numericInput("nMin", label = "Minimum N counts per sample", value = 50),
                                numericInput("prc", label = "Prior count", value = 0.5)
                         ),
                         column(3, 
                                h4("ParameterB"),
                                numericInput("minSpan", label = "Minimum span", value = 4.0)
                         )
                     )
    ),
    conditionalPanel(condition="input.tab == 1",
                     fluidRow(
                         column(9,
                                DT::dataTableOutput('datab')
                                ),
                         column(3, 
                                uiOutput("selDataSet1"),
                                numericInput("nMin1", label = "Minimum N counts per sample", value = 50),
                                numericInput("prc1", label = "Prior count", value = 0.5),
                                numericInput("minSpan1", label = "Minimum span", value = 4.0),
                                uiOutput("selCol"),
                                column(6, verbatimTextOutput('tout'))
                                
                         )
                     )
    )
    )
    
    
)

