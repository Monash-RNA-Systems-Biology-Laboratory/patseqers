library(shiny)
source("helper.R")


shinyUI(fluidPage(
    
    # Application title
    titlePanel("Genewise Plotter"),
    sidebarPanel(
        uiOutput("experi_sel"),  
        
        selectInput("org.type", label = "Type of organism", choices =
                        list("Humans" = "hsapiens_gene_ensembl", "Yeast" = "scerevisiae_gene_ensembl", 
                             "Ce.elegans" = "celegans_gene_ensembl"), 
                    selected = "hsapiens_gene_ensembl"),
        helpText("This is required to ")
        checkboxInput("x.rep", label = "Plot means of selected replicates on the x-axis", 
                      value = TRUE),
        uiOutput("x.sel.ui"),
        selectInput("datatypex", label = "Type of data to be plotted on x-axis", choices =
                        list("Gene expression (RPM)" = 1, "Tail length" = 2, "Peak pair shift" = 3), selected = 1),
        checkboxInput("y.rep", label = "Plot means of selected replicates on the y-axis", 
                      value = TRUE),
        uiOutput("y.sel.ui"),
        selectInput("datatypey", label = "Type of data to be plotted on y-axis", choices =
                        list("Gene expression (RPM)" = 1, "Tail length" = 2,  "Peak pair shift" = 3), selected = 1),
        numericInput("numfilter", label = "Minimum count filter", value = 25, min = 0),
        textInput("file_name", label = "Name of file to download", "file"),
        downloadButton("deps", label = "Download eps file"),
        downloadButton("dpdf", label = "Download pdf file")
        
    ),
    
    
    mainPanel(
        textOutput("total.txt"),
        plotOutput("gplot"),
        radioButtons("type_select", "Choose whether to search for genes of interest or 
                     paste a list of genes (red)", choices= list("Paste genes" = 1, "Search genes (this 
                                                            may take a while to load and will 
                                                                 initially give an error message)" = 2), 
                     selected = 1),
        conditionalPanel(
            condition = "input.type_select == 1",
            textInput("gene.paste.sel", label = "Paste genes to highlight", value = "")
        ),
        conditionalPanel(
            condition = "input.type_select == 2",
            uiOutput("gene.search.ui")
        ),
        textOutput("sel.txt"),
        textInput("key.select", label = "Search a key term (blue)", value = ""),
        textOutput("key.txt"),
        textInput("goterm.select", label = "Search a GOTerm (green)", value = ""),
        textOutput("goterm.txt"),
        tableOutput("gene.table")
        
    )
)
)
