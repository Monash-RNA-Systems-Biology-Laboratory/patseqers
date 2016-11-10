library(shiny)
source("helper.R")


shinyUI(fluidPage(
  
  # Application title
  titlePanel("Poly A tail and gene expression plotter"),
  sidebarPanel(
    radioButtons("radio.x", label = h4("Select condition for the x-axis"),
                 choices = list("N2" = "N2", "Gld2" = "Gld2", "Cpb3" = "Cpb3", "Gld2pcb19" = "Gld2_pcb19", 
                                "Gld2-CCF1" = "Gld2_CCF1", "X1.2.cell.egg" = "X1.2.cell.egg"),
                 selected = "N2", inline = T
    ),
    radioButtons("datatypex", label = "Type of data to be plotted on x-axis", choices =
                                 list("Gene expression (RPM)" = "gene expression",
                                      "Tail length" = "tail length"),
                 selected = "tail length", inline = T
    ),
    radioButtons("radio.y", label = h4("Select condition for the y-axis"),
                 choices = list("N2" = "N2", "Gld2" = "Gld2", "Cpb3" = "Cpb3", "Gld2-pcb19" = "Gld2_pcb19", 
                                "Gld2-CCF1" = "Gld2_CCF1", "X1.2.cell.egg" = "X1.2.cell.egg"),
                 selected = "Gld2", inline = T
    ),
    radioButtons("datatypey", label = "Type of data to be plotted on y-axis", choices =
                   list("Gene expression (RPM)" = "gene expression",
                        "Tail length" = "tail length"),
                 selected = "tail length", inline = T
    ),
    
    # selectInput("x.selection", label = h4("Select samples for the x-axis"),
    #             choices = specific_samples,
    #             selected =  specific_samples[1:3], multiple = TRUE),
    # selectInput("datatypex", label = "Type of data to be plotted on x-axis", choices =
    #               list("Gene expression (RPM)" = "gene expression",
    #                    "Tail length" = "tail length"),
    #             selected = "tail length"),
    # selectInput("y.selection", label = h4("Select samples for the y-axis"),
    #             choices = specific_samples,
    #             selected =  specific_samples[4:6], multiple = TRUE),
    # selectInput("datatypey", label = "Type of data to be plotted on y-axis", choices =
    #               list("Gene expression (RPM)" = "gene expression", 
    #                    "Tail length" = "tail length"), 
    #             selected = "tail length"),
    numericInput("numfilter", label = "Minimum count filter", value = 100, min = 0),
    radioButtons("type_filter", "Choose whether to filter by just one selection or all selections for
                     the selected axis", 
                 choices= list("One selection" = 1, "All selections" = 2), 
                 selected = 1),
    helpText("The filter will be applied to both datatypes. It applies a filter that requires
                 at least one of the selected samples to meet the filter or all selected samples meets the set filter for a 
                 gene to be retained in the dataset. Selected samples refers to the total samples across the X and Y axis."),
    textInput("file_name", label = "Name of graph to download", "file"),
    downloadButton("deps", label = "Download eps file"),
    downloadButton("dpdf", label = "Download pdf file")
  ),
  
  
  mainPanel(
    tabsetPanel(type = "tabs",
                tabPanel("Plot", 
                         ##sidebarLayout(
                          # mainPanel(
                         # div(
                         #   style = "position:relative",
                         plotOutput("gplot", width = "550px", height = "550px", 
                                    # hover=hoverOpts(id="p_hover", delay = 100, delayType = "debounce")),
                                    click = "plot_click"),
                         # uiOutput("hover_over_info"),
                         # width = 7
                         #),
                         helpText(""),
                         textOutput("total.txt"),
                         verbatimTextOutput("info_plot"),
                         textOutput("sel.txt"),
                         textOutput("key.txt"),
                         textOutput("goterm.txt"),
                         tableOutput("debug.goterm")
                         #  ),
                         #sidebarPanel(position = ("right"),
                        #   helpText("Testing sidebar panel")
                        # ) 
                     #    )
                ),
                tabPanel("Search genes",
                         helpText(h4("After selecting the genes and GO Terms of interest, return to the plot tab and the plot will update itself.
                                     If you would like to download the searched genes, go to the Output Tables tab.")),
                         radioButtons("type_select", "Choose whether to search for genes of interest or 
                     paste a list of genes (purple)", choices= list("Paste a list of genes" = 1, "Find genes" = 2), 
                                      selected = 1),
                         conditionalPanel(
                           condition = "input.type_select == 1",
                           helpText("Note: Seperate each gene by a space. Pasting a list of genes requires the genes to match exactly the 
                                names on the file. The dataset uses RefSeq identifiers (e.g NM_058260), gene (e.g homt-1) 
                                    and transcript names (e.g Y74C9A.3). Every gene should have a RefSeq ID but varies on whether the gene name
                                    or transcript name was assigned to it. The link below is to a tabs file that contains all the genes in the dataset."),
                           tags$a(href='data/gld2_dataset_names.txt', target='blank', 'List of all the genes in this dataset', download = 'gld2_dataset_names.txt'),
                           textInput("gene.paste.sel", label = "Paste genes to highlight", 
                                     value = "mex-5 mex-6 puf-3 cbd-1 pos-1 mex-3 oma-1 air-1 C05C10.5 puf-7 puf-11")
                         ),
                         conditionalPanel(
                           condition = "input.type_select == 2",
                           helpText("Only 10 gene names will be displayed at a time but you can type into the 
                                          box to find the gene you are interested it. If you can't find the gene you are looking for, it may not be in the dataset or the filter may need to be adjusted"),
                           uiOutput("gene.search.ui")
                         ),
                         
                         
                         radioButtons("type_goterm", "Choose whether to paste a GoTerm or find a GO slim accession",
                                      choices = list("Paste GOTerm" = 1, "Search GO slim" = 2), selected = 2),
                         conditionalPanel(
                           condition = "input.type_goterm == 1",
                           textInput("goterm.select", label = "Paste a GOTerm (blue)", value = ""),
                           helpText("Only input one GO term.")
                         ),
                         conditionalPanel(
                           condition = "input.type_goterm == 2",
                           selectizeInput("goslim.select", label = "Select GOslim accession",
                                          choices = goslim_list,
                                          selected =  "GO:0003735", multiple = FALSE)
                         ),
                         helpText("If the following error appears: 'Error 1: Extra content at the end of the document' on the Plot tab then the
                                      Biomart servers are down and the GO Term/Slim search will not be useable")
                ),
                tabPanel("Output Tables",
                         helpText("Selected genes"),
                         tableOutput("sel.table"),
                         downloadButton("dsel", label = "Download selected genes table"),
                         helpText("GO term genes"),
                         tableOutput("goterm.table"),
                         downloadButton("dgot", label = "Download table for GO term genes"),
                         tableOutput("annotated_selected")
                )
    )
  )
)
)
