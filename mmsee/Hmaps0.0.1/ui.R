#ui.R
#Authors - Michael See    mmsee2@student.monash.edu.au

shinyUI(fluidPage(
    
    title = "Heatmap",
    tabsetPanel(
        tabPanel("Options",
                 #     plotOutput("plot1", width="1400", height="850"), #Testing the output
                 fluidRow(
                     column(3,
                            radioButtons("tab", label="Use Table", choices=list("False"=0, "True"=1), selected = 0, inline = T),
                            h6("Please note that libraries may take up to a minute to load on upon launch.")
                     ),
                     column(3,
                              numericInput("dwidth", label = "Plot Width (in)", value = 14)
                     ),
                     column(3,
                            numericInput("dheight", label = "Plot Height (in)", value = 20)
                     ),
                     column(3,
                            downloadButton("download", "Download the current plot")
                     )),
                 
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
                                             h4("Plot Information"),
                                             numericInput("minSpan", label = "Minimum span", value = 4.0),
                                             uiOutput("selCol")
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
                                             uiOutput("selCol1"),
                                             column(6, verbatimTextOutput('tout'))
                                             
                                      )
                                  )
                 )
        ),
        tabPanel("Heatmap", 
                 h2("Heatmap"),
                 hr(),
                 uiOutput("plotui"),
                 fluidRow(
                     column(3,
                            numericInput("pwidth", label = "Plot Width (px)", value = 1400)
                     ),
                     column(3,
                            numericInput("pheight", label = "Plot Height (px)", value = 800)
                     ),
                     column(3,
                            actionButton("gopt", "Resize Plot")
                            #verbatimTextOutput("brushOut")
                     )
                     )
                     
                 
        
    )
)

)
)
