source("major_calculations.R")
library("ggplot2")
library(reshape)
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
  groupings <-reactive({
    lapply(find_bam_files(bam_file_path), function(i) input[[paste0('snumber', i)]])    
  })
  
  poly_a_counts<- reactive({
 
    initial_table <- get_a_counts (bam_file_path, gffInput(),input$select_bam_files,
                                   input$select_genes,groupings())
    subsetted_by_sliders <- initial_table[initial_table$number_of_ad_bases >= 
                                            input$ad_slider&
                                            initial_table$width >=
                                            input$al_length[1]&
                                            initial_table$width <=
                                            input$al_length[2],]
  })
  output$print_poly_a_counts <- renderDataTable({
    poly_a_counts()    
  })
  output$n_reps<- renderText ({
    input$n_replicates   
  })

  output$gene_info <- renderText({
    full_df <- poly_a_counts()  
    split_frame <- split(full_df, full_df$sample)
    names_string (split_frame, input$merge) 
  })
  plot_calcs <- reactive({
    make_plot(poly_a_counts(), input$xslider,input$select_genes, input$legend, input$merge)
  })    
  output$scp_plot<- renderPlot({  
      plot_calcs()
  })
  #Workaround for a shiny bug thatdoesn't handle reactive plots well. 
  plot_calcs2 <- function(){
    make_plot(poly_a_counts(), input$xslider,input$select_genes, input$legend)
  }
  output$downloadPlot <- downloadHandler(
    filename = function(){
      paste(trim(input$select_genes), '.pdf', sep='')
    },
    content = function(file){
      pdf(file,width = 10)
      postscript(file)
      plot_calcs2()     
      dev.off()
  
  })
    
  output$selected_dataset<- renderText({
    selected_data(input$dataset)
  })
  
})