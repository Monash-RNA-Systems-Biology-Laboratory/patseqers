library(shiny)
source("helper.R")

shinyUI(pageWithSidebar(
  headerPanel(h3("Changes in gene expression after C.albicans infection")),
  
  sidebarPanel(
      selectInput("select.experiment", label = "Select experiment", 
                  choices = list.dirs(full.names=F, recursive =F), 
                  selected =  list.dirs(full.names=F, recursive =F)[2]),
      radioButtons("org.sel", "Choose which organism to plot", 
                   choices= list("C.albicans" = 1, "Mouse" = 2), 
                   selected = 1),
      helpText("This is important when changing the dataset to look at. If one of the samples is not 
               displaying on the graph, check that the correct organism is selected"),
      numericInput("num.filter", label = "Minimum count filter", value = 25, min = 0),
      uiOutput("gene.select"),
      radioButtons("plot.view", "View just the gene or the gene of interest in comparison with expression changes
                   overall", 
                   choices= list("Only gene of interest" = 1, "Comparison with overall expression change" = 2), 
                   selected = 1),
      downloadButton("download.eps", label = "Download eps file")
  ),
  
  mainPanel(
     plotOutput("time.plot", width = "600px", height = "600px", click = "plot_click"),
     #tableOutput("test.table"),
     verbatimTextOutput("info.txt"),
     helpText("You can click on the graph to get X & Y values. Note that the only measured data occur at the given
              timepoints."),
     helpText("Gene expressions refers to the mean of the normalised counts for the replicates of either condition 
              at each time point. The varistran package was used to normalised the raw counts.")
    
  )
  
  
))
