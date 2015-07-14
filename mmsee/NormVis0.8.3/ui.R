#ui.R

shinyUI(fluidPage(
    titlePanel('Visualisation of Normalised vs Un-Normalised Data'),
    sidebarLayout(
        sidebarPanel(
            uiOutput("selDataSet"),
            
            
            numericInput("nMin", label = "Minimum N counts per sample", value = 50),
            uiOutput("SelGeneNorm"), #Add back uiOutput("SelgeneGraph")
            
            
            numericInput("pRC", label = "Prior Count for Normalisation against Library", value = 0.5),
            uiOutput("selGene"),
	    checkboxInput("drTMM", "Draw TMM instead of naive normalisation", FALSE),
            actionButton("do", "Enable interactivity")
            
        ),
        mainPanel(
            
            plotOutput('plot2'),
            plotOutput('plot3'),
	    downloadButton('PlotD1', 'Download Plot 1'),
	    downloadButton('PlotD2', 'Download Plot 2')
		            
            
        )
    )
))
