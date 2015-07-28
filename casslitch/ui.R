library(shiny)

shinyUI(fluidPage(
  
  # Application title
  titlePanel("Genewise Plotter"),
  
  sidebarPanel(
    selectInput("sampleX", "Select sample for x-axis",
                       check_boxes,selected="N2_mean"),
    
    selectInput("sampleY", "Select sample for y-axis",
                       check_boxes,selected="Gld2_mean"),
    
    radioButtons("dfxtype", "Data type x-axis", 
                       choices= list("expression log2(RPM)"= "genewise_exp", "tail length"="genewise_tail_length"),selected="tail length"),
    
    radioButtons("dfytype", "Data type y-axis", 
                       choices= list("expression log2(RPM)"= "genewise_exp", "tail length"="genewise_tail_length"),selected="tail length"),
    
    radioButtons("isoforms", "Include isoforms in GO term search?", 
                       choices= list("Yes"= T, "No"= F),selected="No")
  ),
  
  
  mainPanel(
      plotOutput("distPlot"),
      submitButton(text = "Apply Changes", icon = NULL),
      textInput("bygene", label="search by gene name", value = ""),
      textInput("GOterm", label="search by GOterm", value = ""),
      textInput("productterm", label="search by product description key term", value = ""),
      downloadButton("downloadPlot", label = "Download Plot")
  
    )
  )
)
