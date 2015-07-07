source("major_calculations.R")
library(Rsamtools)

shinyServer(function(input, output, session) {
  
  # Combine the selected variables into a new data frame
  selectedData <- reactive({
    input$dataset
  })
  number_of_datasets <- reactive({
    input$dataset
  })  
  read_files <- reactive({
    read_required_files(input$dataset)
  })
  dataset1 <- reactive({
    read_required_files(input$dataset)
  }) 
  gffInput <- reactive({
  filter_gff_for_rows(paste(gff_file_path,input$gff_select,sep=""), input$select_genes)
  })
  output$gff_rows<- renderDataTable({
    gffInput () 
  })
  poly_a_counts<- reactive({
    get_a_counts (bam_file_path, gffInput(),input$select_bam_files,input$select_genes)
    
  })
  output$print_poly_a_counts <- renderDataTable({
    poly_a_counts() 
    
  })
  output$scp_plot<- renderPlot({
    make_plot()
    
  })
  output$selected_dataset<- renderText({
    selected_data(input$dataset)
  })
  
})