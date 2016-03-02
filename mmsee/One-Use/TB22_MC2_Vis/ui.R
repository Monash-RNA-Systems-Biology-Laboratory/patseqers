
shinyUI(fluidPage(
    headerPanel("Plot1_incomplete"),
    sidebarPanel(width=3,
                 uiOutput('selgene'),
                 uiOutput('normdisp'),
                 downloadButton('PlotD1', 'Download Plot')
                 ),
    mainPanel(plotOutput('plot'))
    ))