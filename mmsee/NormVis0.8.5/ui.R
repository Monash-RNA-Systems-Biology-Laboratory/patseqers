#ui.R
#Authors: Michael See

shinyUI(fluidPage(
    titlePanel('Visualisation of Normalised data'),
    sidebarLayout(
        sidebarPanel(
            uiOutput("selDataSet"),
            
            
            numericInput("nMin", label = "Minimum N counts per sample", value = 50),
            uiOutput("SelGeneNorm"), 
            
            
            numericInput("pRC", label = "Prior Count for Normalisation against Library", value = 0.5),
            uiOutput("selGene"),
	        checkboxInput("drTMM", "Draw TMM instead of normalisation against 1 gene", FALSE),
            actionButton("do", "Enable interactivity"),
            verbatimTextOutput('tout')
        
            
        ),
        mainPanel(
            
            plotOutput('plot2'),
            plotOutput('plot3'),
	    downloadButton('PlotD1', 'Download Plot 1'),
	    downloadButton('PlotD2', 'Download Plot 2')
		            
            
        )
    )
))
