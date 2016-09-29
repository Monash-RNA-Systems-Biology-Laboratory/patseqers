library(shiny)
source("helper.R")


shinyUI(fluidPage(
  
  # Application title
  titlePanel("Poly A tail and gene expression plotter"),
  sidebarPanel(
    selectInput("x.selection", label = h4("Select samples for the x-axis"),
                choices = specific_samples,
                selected =  specific_samples[1:3], multiple = TRUE),
    selectInput("datatypex", label = "Type of data to be plotted on x-axis", choices =
                  list("Gene expression (RPM)" = "gene expression", 
                       "Tail length" = "tail length"), 
                selected = "Gene expression"),
    selectInput("y.selection", label = h4("Select samples for the y-axis"),
                choices = specific_samples,
                selected =  specific_samples[4:6], multiple = TRUE),
    selectInput("datatypey", label = "Type of data to be plotted on y-axis", choices =
                  list("Gene expression (RPM)" = "gene expression", 
                       "Tail length" = "tail length"), 
                selected = "Gene expression"),
    numericInput("numfilter", label = "Minimum count filter", value = 100, min = 0),
    radioButtons("type_filter", "Choose whether to filter by just one selection or all selections for
                     the selected axis", 
                 choices= list("One selection" = 1, "All selections" = 2), 
                 selected = 1),
    helpText("The filter will be applied to both datatypes. It applies a filter that requires
                 at least one sample to meet the filter or all samples meets the set filter for a 
                 gene to be retained in the dataset."),
    textInput("file_name", label = "Name of file to download", "file"),
    downloadButton("deps", label = "Download eps file"),
    downloadButton("dpdf", label = "Download pdf file")
  ),
  
  
  mainPanel(
    tabsetPanel(type = "tabs",
                tabPanel("Plot", 
                         
                         plotOutput("gplot", width = "550px", height = "550px", click = "plot_click"),
                         helpText(""),
                         textOutput("total.txt"),
                         verbatimTextOutput("info_plot"),
                         textOutput("sel.txt"),
                         textOutput("key.txt"),
                         textOutput("goterm.txt"),
                         tableOutput("debug.goterm")
                ),
                tabPanel("Data select",
                         helpText(h4("After selecting the genes and GO Terms of interest, return to the plot tab and the plot will update itself")),
                         radioButtons("type_select", "Choose whether to search for genes of interest or 
                     paste a list of genes (purple)", choices= list("Paste a list of genes" = 1, "Find genes" = 2), 
                                      selected = 1),
                         conditionalPanel(
                           condition = "input.type_select == 1",
                           helpText("Note:Seperate each gene by a space. Pasting a list of genes requires the genes to match exactly the 
                                names on the file. The dataset uses RefSeq identifies (e.g NM_058260), gene (e.g homt-1) 
                                    and transcript names (e.g Y74C9A.3). Every gene should have a RefSeq ID but varies on whether the gene name
                                    or transcript name was assigned to it."),
                           tags$a(href='data/gld2_dataset_names.txt', target='blank', 'List of all the genes in this dataset', download = 'gld2_dataset_names.txt'),
                           textInput("gene.paste.sel", label = "Paste genes to highlight", value = "homt-1 rcor-1 rab-11.1")
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
                                          selected =  NULL, multiple = FALSE)
                         ),
                         helpText("If the following error appears: 'Error 1: Extra content at the end of the document' on the Plot tab then the
                                      Biomart servers are down and the GO Term/Slim search will not be useable")
                ),
                tabPanel("Output tables",
                         helpText("Selected genes"),
                         tableOutput("sel.table"),
                         downloadButton("dsel", label = "Download selected genes table"),
                         helpText("GO term genes"),
                         tableOutput("goterm.table"),
                         downloadButton("dgot", label = "Download table for GO term genes")
                )
    )
  )
)
)