library(shiny)

shinyUI(fluidPage(
  titlePanel("Gene cross infection app"),
  
  sidebarLayout(
    sidebarPanel(
      # selectInput("select_experiment", label = h3("Select experiment"),
      #             choices = list.dirs(full.names=F, recursive =F),
      #             selected =  list.dirs(full.names=F, recursive =F)[1]),
      numericInput("num_filter", label = "Minimum count filter", value = 25, min = 0),
      uiOutput("gene_select"),
      radioButtons("plot_mode", "Select which mode to view",
                   choices= list("Only gene of interest" = 1,
                                 "Comparison with overall expression change" = 2,
                                 "View YPD values for gene of interest" = 3),
                   selected = 1),
      helpText("The plotted line graphs only take into account the time-points from 2.00 to 6.00. To see the YPD values (time 0), select the third option.
               To keep the x-axis limits consistent with the other two view modes, the YPD values are plotted horizontally. 
               E1 refers to MC_119 and E2 refers to samples from MC2_57_59."),
      downloadButton("download_eps", label = "Download eps file")
    ),
    mainPanel (
      tabsetPanel(type = "tabs",
                  tabPanel("Plot", 
                           plotOutput("time_plot", width = "950px", height = "450px", click = "plot_click"),
                           verbatimTextOutput("info_txt"),
                           helpText("You can click on the graph to get X & Y values. Note that the only measured data occur at the given
              timepoints."),
                           helpText("Gene expressions refers to the mean of the normalised counts for the samples of either condition 
              at each time point. The varistran package was used to normalised the raw counts")
                  ),
                  tabPanel("Table of data plotted",
                           tableOutput("table_plotted"),
                           helpText("Experiment 1 refers to MC_119 and experiment 2 refers to samples from MC2_57_59.
                                    Time points between 2 and 3 have been modified in order to be plotted neatly on the x-axis.")
                  ),
                  tabPanel("Library size",
                           tableOutput("library_size_table"),
                           helpText("Library size will vary depending on the count filter")
                           )
      )
      
    ))
))
