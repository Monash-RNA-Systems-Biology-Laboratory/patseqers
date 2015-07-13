library(shiny)

# gene_exp<-read.csv("genewise-count.csv")
# check_boxes<-names(gene_exp[c(2:length(names(gene_exp)))])


# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("TB22-88"),
  
  sidebarPanel(
    selectInput("sampleX", "xvar",
                       check_boxes,selected="WT.rep1"),
    
    selectInput("sampleY", "yvar",
                       check_boxes,selected="WT.rep2")
  ),
  
  
  mainPanel(
      plotOutput("distPlot"),
      textInput("bygene", label="search custom", value = "YAL067C"),
      textInput("GOterm", label="GOterm", value = "GO:0005634")
    )
  )
)
