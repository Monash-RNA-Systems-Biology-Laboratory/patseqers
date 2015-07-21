library(shiny)

# gene_exp<-read.csv("genewise-count.csv")
# check_boxes<-names(gene_exp[c(2:length(names(gene_exp)))])


# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("TB22-88"),
  
  sidebarPanel(
    selectInput("sampleX", "xvar",
                       check_boxes,selected="N2_mean"),
    
    selectInput("sampleY", "yvar",
                       check_boxes,selected="Gld2_mean"),
    
    radioButtons("dataframx", "Data type X", 
                       choices= list("expression log2(RPM)"= "genewise_exp", "tail length"="genewise_tail_length"),selected="tail length"),
    
    radioButtons("dataframy", "Data type Y", 
                       choices= list("expression log2(RPM)"= "genewise_exp", "tail length"="genewise_tail_length"),selected="tail length"),
    radioButtons("isoforms", "Include isoforms in GO term search?", 
                       choices= list("Yes"= T, "No"= F),selected="No")
  ),
  
  
  mainPanel(
      plotOutput("distPlot"),
      submitButton(text = "Apply Changes", icon = NULL),
      textInput("bygene", label="search by gene", value = ""),
      textInput("GOterm", label="GOterm", value = ""),
      textInput("productterm", label="Product description key term", value = ""),
      downloadButton("downloadPlot", label = "Download Plot")
  
    )
  )
)
