library(shiny)
source("helper.R")


shinyUI(fluidPage(
    
    # Application title
    titlePanel("Genewise Plotter"),
    sidebarPanel(
#         uiOutput("experi_sel"),  
        selectInput("select.experiment", label = "Select experiment", 
                    choices = list.dirs(full.names=F, recursive =F), 
                    selected =  list.dirs(full.names=F, recursive =F)[1]),
        selectInput("org.type", label = "Type of organism", choices =
                        list("Human" = "hsapiens_gene_ensembl", "Yeast" = "scerevisiae_gene_ensembl", 
                             "Ce.elegans" = "celegans_gene_ensembl", "Mouse" = "mmusculus_gene_ensembl"), 
                    selected = "celegans_gene_ensembl"),
        helpText("The correct organism only needs to be selected for the GO Term search to work"),
        uiOutput("x.sel.ui"),
        selectInput("datatypex", label = "Type of data to be plotted on x-axis", choices =
                        list("Gene expression (RPM)" = "gene expression", 
                             "Tail length" = "tail length", "Peak pair shift" = "peak pair shift"), 
                    selected = "Gene expression"),
        uiOutput("y.sel.ui"),
        selectInput("datatypey", label = "Type of data to be plotted on y-axis", choices =
                        list("Gene expression (RPM)" = "gene expression", 
                             "Tail length" = "tail length", "Peak pair shift" = "peak pair shift"), 
                    selected = "Gene expression"),
        numericInput("numfilter", label = "Minimum count filter", value = 25, min = 0),
        radioButtons("type_filter", "Choose whether to filter by one selection per axis or all selections for
                     the selected axis", 
                     choices= list("One selection" = 1, "All selections" = 2), 
                     selected = 1),
        helpText("The filter will be applied to all three datatypes. It applies a filter that requires
                 at least one sample out of the chosen samples for each axis or all samples meets the set filter for a 
                 gene to be retained in the dataset."),
        textInput("file_name", label = "Name of file to download", "file"),
        downloadButton("deps", label = "Download eps file"),
        downloadButton("dpdf", label = "Download pdf file"),
        radioButtons("plot_axes", "Choose to standardise the axes range", choices= 
                         list("Unstardardise (axes range might be unequal)" = 1, "Standardise (identical range) - not
                              recommended for different data types" = 2), 
                     selected = 1)
        
    ),
    
    
    mainPanel(
        tabsetPanel(type = "tabs",
                    tabPanel("Plot", 
                             
                             plotOutput("gplot", width = "550px", height = "550px", click = "plot_click"),
                             textOutput("total.txt"),
                             verbatimTextOutput("info_plot"),
                             textOutput("sel.txt"),
                             textOutput("key.txt"),
                             textOutput("goterm.txt")
                    ),
                    tabPanel("Data Select",
                             radioButtons("type_select", "Choose whether to search for genes of interest or 
                     paste a list of genes (purple)", choices= list("Paste genes" = 1, "Search genes" = 2), 
                                          selected = 1),
                             helpText("Note: Pasting a list of genes requires the genes to match exactly the 
                                names on the file. Seperate each gene by a space"),
                             conditionalPanel(
                                 condition = "input.type_select == 1",
                                 textInput("gene.paste.sel", label = "Paste genes to highlight", value = "")
                             ),
                             conditionalPanel(
                                 condition = "input.type_select == 2",
                                 helpText("It may take the search input a while to load"),
                                 uiOutput("gene.search.ui")
                             ),
                             
                             textInput("key.select", label = "Search a key term (yellow)", value = ""),
                             
                             helpText("Only input one key term"),
                             textInput("goterm.select", label = "Search a GOTerm (blue). Remember to
                                       check which organism has been selected", value = ""),
                             
                             helpText("Only input one GO term")
                    ),
                    tabPanel("Selection table",
                             helpText("Selected genes"),
                             tableOutput("sel.table"),
                             downloadButton("dsel", label = "Download selected genes table"),
                             helpText("Key term genes"),
                             tableOutput("key.table"),
                             downloadButton("dkey", label = "Download table for key term genes"),
                             helpText("GO term genes"),
                             tableOutput("goterm.table"),
                             downloadButton("dgot", label = "Download table for GO term genes")
                             )
        )
    )
)
)
