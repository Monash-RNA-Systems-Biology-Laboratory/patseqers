library(shiny)

shinyUI(fluidPage(
  
  # Application title
  titlePanel("Genewise Plotter"),
  
  sidebarPanel(
    uiOutput("select_file_path"),
    selectInput("sampleX", "Select sample for x-axis",
                       check_boxes,selected="N2_mean"),
    
    selectInput("sampleY", "Select sample for y-axis",
                       check_boxes,selected="Gld2_mean"),
    
    radioButtons("dfxtype", "Data type x-axis", 
                       choices= list("expression log2(RPM)"= "genewise_exp", "tail length"="genewise_tail_length"),selected="tail length"),
    
    radioButtons("dfytype", "Data type y-axis", 
                       choices= list("expression log2(RPM)"= "genewise_exp", "tail length"="genewise_tail_length"),selected="tail length"),
    
    radioButtons("isoforms", "Include isoforms in GO term search?", 
                       choices= list("Yes"= T, "No"= F),selected="No"),
    radioButtons("pval", "Show p value on chart?", 
                 choices= list("Yes"= T, "No"= F),selected="Yes"),
    radioButtons("organism", "Select an organism", 
                 choices= list("Worm"= "celegans_gene_ensembl", "Yeast"= "scerevisiae_gene_ensembl"),selected="celegans_gene_ensembl")

  ),
  
  
  mainPanel(
      plotOutput("distPlot"),
      textOutput("print"),
      submitButton(text = "Apply Changes", icon = NULL),
      textInput("bygene", label="search by WBID (seperate IDs by a single space)", value = ""),
      textInput("byyeastgene", label="search by yeast gene name (or start of gene name)", value = ""),
      textInput("GOterm", label="search by GOterm (seperate terms by ,)", value = ""),
      textInput("productterm", label="search by product description key term (seperate terms by ,)", value = ""),
      downloadButton("downloadPlot", label = "Download Plot")
  
    )
  )
)
