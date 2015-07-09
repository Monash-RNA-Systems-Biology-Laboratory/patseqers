#ui.R

shinyUI(fluidPage(
    titlePanel('Visualisation of Normalised vs Un-Normalised Data'),
    sidebarLayout(
        sidebarPanel(
            
            numericInput("nMin", label = "Minimum N counts per sample", value = 50),

            uiOutput("SelGeneNorm"), 

            uiOutput("selGene"),
            
            numericInput("pRC", label = "Prior Count for Normalisation against Library", value = 0.5),
            
            actionButton("do", "Submit and Plot")
            
            
        ),
        mainPanel(
            
            plotOutput('plot2'),
            plotOutput('plot3')

            

        )
    )
))