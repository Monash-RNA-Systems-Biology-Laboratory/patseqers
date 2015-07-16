library(shiny)


shinyUI(fluidPage(
  
  # Application title
  titlePanel("TB22-88"),
  
  sidebarPanel(
    selectInput("sampleX", "x-axis",
                       check_boxes,selected="WT.rep1"),
    
    selectInput("sampleY", "y-axis",
                       check_boxes,selected="WT.rep2"),
    
    radioButtons("dataframx", "Data type X-axis", 
                       choices= list("expression log2(RPM)"= "genewise_exp", "tail length"="genewise_tail_length"),selected="expression log2(RPM)"),
    
    radioButtons("dataframy", "Data type Y-axis", 
                       choices= list("expression log2(RPM)"= "genewise_exp", "tail length"="genewise_tail_length"),selected="expression log2(RPM)"),
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
