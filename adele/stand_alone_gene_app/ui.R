library(shiny)
library(shinyBS)
source("helper.R")
load("www/data/preload_goslim_list.rds")

shinyUI(fluidPage(
  
  # Application title
  titlePanel("Poly A tail and gene expression plotter"),
  sidebarPanel(
    h4("Options"),
    selectInput("input_x", label = "Select condition for the x-axis",
                choices = list("N2" = "N2", "gld2" = "gld2", "gld2-pcb19-RNAi" = "gld2_pcb19_RNAi", 
                               "gld2-ccf1-RNAi" = "gld2_ccf1_RNAi", "1-2-cell" = "1.2_cell"),
                selected = "N2", width = "80%"
    ),
    bsTooltip(id = "input_x", title = "Select which strain to plot on the x-axis ", 
              placement = "top", trigger = "hover"),
    radioButtons("datatypex", label = "Type of data to be plotted on x-axis", choices =
                   list("Gene expression (RPM)" = "gene expression",
                        "Tail length" = "tail length"),
                 selected = "tail length", inline = T
    ),
    bsTooltip(id = "datatypex", title = "Select whether to plot gene expression (RPM) or tai length data on the x-axis", 
              placement = "top", trigger = "hover"),
    selectInput("input_y", label = "Select condition for the y-axis",
                choices = list("N2" = "N2", "gld2" = "gld2", "gld2-pcb19-RNAi" = "gld2_pcb19_RNAi", 
                               "gld2-ccf1-RNAi" = "gld2_ccf1_RNAi", "1-2-cell" = "1.2_cell"),
                selected = "gld2", width = "80%"
    ),
    bsTooltip(id = "input_y", title = "Select which strain to plot on the y-axis ", placement = "top", trigger = "hover"),
    radioButtons("datatypey", label = "Type of data to be plotted on y-axis", choices =
                   list("Gene expression (RPM)" = "gene expression",
                        "Tail length" = "tail length"),
                 selected = "tail length", inline = T
    ),
    bsTooltip(id = "datatypey", title = "Select whether to plot gene expression (RPM) or tai length data on the y-axis",
              placement = "top", trigger = "hover"),
    
    tags$hr(),
    
    checkboxInput(inputId = "advanced_options", label = "Display advanced options for filtering", value = FALSE),
    conditionalPanel(
      condition = "input.advanced_options == true",
      numericInput("numfilter", label = "Minimum count filter", value = 100, min = 0),
      bsTooltip(id = "numfilter", 
                title = "Select the minimum number of counts for filtering",
                placement = "top", trigger = "hover"),
      radioButtons("type_filter", "Choose whether to filter by just one sample or all sample for the selected conditions", 
                   choices= list("One sample" = 1, "All samples" = 2), 
                   selected = 1, inline = T),
      bsTooltip(id = "type_filter", title = "Help: The filter will be applied to both datatypes. For one selection - across the two selected conditions, at least one sample must meet the filter. For all samples - across the two selected conditions, all samples must meet the filter.",
                placement = "top", trigger = "hover")
    ),
    
    tags$hr(),
    
    h4("Search Genes"),
    
    radioButtons("type_goterm", "Choose whether to select a GO slim accession or paste a GoTerm",
                 choices = list("Select GO slim" = 1, "Paste GOTerm" = 2), selected = 1, inline = T),
    bsTooltip(id = "type_goterm", title = 
                "If the following error appears: 'Error 1: Extra content at the end of the document' on the Plot tab then the Biomart servers are down and the GO Term/Slim search will not be useable",
              placement = "top", trigger = "hover"),
    conditionalPanel(
      condition = "input.type_goterm == 1",
      selectizeInput("goslim_select", label = "Find GOslim accession",
                     choices = slim_list,
                     selected =  "GO:0003735", multiple = FALSE, width = "80%"),
      bsTooltip(id = "goslim_select", title = "GO Slims are cut-down subset of Gene Ontologies",
                placement = "top", trigger = "hover")
    ),
    conditionalPanel(
      condition = "input.type_goterm == 2",
      textInput("goterm_select", label = "Input a GOTerm", value = "", width = "80%"),
      bsTooltip(id = "goterm_select", title = "Type in only one GOTerm. Searching by GOTerm can be more specific than using GO Slims.", 
                placement = "top", trigger = "hover")
    ),
    
    tags$br(),
    
    radioButtons("type_select", "Choose whether to select genes of interest or paste a list of genes", choices= list("Find genes" = 1, "Paste a list of genes" = 2), 
                 selected = 2, inline = T),
    conditionalPanel(
      condition = "input.type_select == 1",
      uiOutput("gene_search_ui")
      
    ),
    conditionalPanel(
      condition = "input.type_select == 2",
      
      textInput("gene_paste_sel", label = "Input genes to highlight", 
                value = "mex-5 mex-6 puf-3 cbd-1 pos-1 mex-3 oma-1 air-1 C05C10.5 puf-7 puf-11", width = "80%"),
      bsTooltip(id = "gene_paste_sel", title = "Help: Seperate each gene by a space. Pasting a list of genes requires the genes to match exactly the names on the file. The dataset uses RefSeq identifiers (e.g NM_058260), gene (e.g homt-1) and transcript names (e.g Y74C9A.3). Every gene should have a RefSeq ID but varies on whether the gene name or transcript name was assigned to it. The link below is to a file that contains all the genes in the dataset.", 
                placement = "top", trigger = "hover"),
      tags$a(href='data/gld2_dataset_names.txt', target='blank', 'List of all the genes in this dataset', download = 'gld2_dataset_names.txt')    )
    
  ),
  
  
  mainPanel(
    tabsetPanel(type = "tabs",
                tabPanel("Plot", 
                         
                         plotOutput("gplot", width = "650px", height = "500px", 
                                    click = "plot_click"),
                         tags$hr(),
                         h4("Plot information"),
                         tableOutput("report_stats"),
                         bsTooltip(id = "report_stats", title = "Information on the number of genes currently being plotted.",
                                   placement = "top", trigger = "hover"),
                         verbatimTextOutput("info_plot"),
                         bsTooltip(id = "info_plot", title = 
                                     "The top 5 closest genes to the clicked mouse position, arranged in order of nearest distance to furthest.",
                                   placement = "top", trigger = "hover"),
                         textInput("file_name", label = "Name of graph to download", "file", width = "50%"),
                         downloadButton("deps", label = "Download eps file"),
                         downloadButton("dpdf", label = "Download pdf file")
                         
                ),
                
                tabPanel("Output Tables",
                         fluidRow(
                           column(4,
                                  helpText("Individual selected genes"),
                                  downloadButton("dsel", label = "Download selected genes"),
                                  tableOutput("sel.table")
                           ),
                           column(4,
                                  helpText("GO Term / Slim annotated genes"),
                                  downloadButton("dgot", label = "Download GO Term/Slim annotated genes"),
                                  tableOutput("annotated_selected")#,
                           )
                         )     
                )
    )
  )
)
)

