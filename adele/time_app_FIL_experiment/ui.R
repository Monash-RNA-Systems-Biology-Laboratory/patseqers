library(shiny)

shinyUI(fluidPage(
  
  # Create a title panel for the app
  titlePanel("mDivi time series app"),
  
  # Use a sidebar layout
  sidebarLayout(
    # Create a sidebar panel for app controls
    sidebarPanel(
      # Create a numeric input to get values for a count filter
      numericInput("num_filter", label = "Minimum count filter", value = 25, min = 0),
      
      # Create a uiOutput to select a gene to plot. This is generated on the server side as the count filter can remove genes with low numbers of genes
      uiOutput("gene_select"),
      
      # Change whether to view the normalised counts of the gene by themselves or view the overall trend
      radioButtons("plot_mode", "Select which mode to view",
                   choices= list("Only gene of interest" = 1,
                                 "Comparison with overall expression change" = 2),
                   selected = 1),
      # Download the current plot as an EPS file
      downloadButton("download_eps", label = "Download eps file")
    ),
    
    # Main panel for containing the plot
    mainPanel (
      
      # Create a plot output for the main plot. Enable the plot to be clicked to provide timepoint and expression information
      plotOutput("time_plot", width = "95%", height = "450px", click = "plot_click"),
      
      # Create an ouput for the click information to be displayed
      verbatimTextOutput("info_txt"),
      
      # Create help text to explain app functionality
      helpText("You can click on the graph to get X & Y values. Note that the only measured data occur at the given
              timepoints."),
      helpText("Gene expression refers to the mean of the normalised log2 CPM for the samples of either condition at each time point. The varistran package was used to normalised the raw counts")
      
    ))
))
