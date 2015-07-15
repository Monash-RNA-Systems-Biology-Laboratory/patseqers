source("major_calculations.R")
library("ggplot2")
library(reshape)
library(Rsamtools)



shinyServer(function(input, output, session) {

  
  # Combine the selected variables into a new data frame
  found_gff_files <- reactive({
    find_gff_files(gff_file_path)
  })
  output$gff_files <- renderUI({
    selectInput("gff_select", label = h4("GFF File Selection"), 
                choices = found_gff_files(), 
                selected =  found_gff_files()[1])   
  })
  found_bam_files <- reactive({
    find_bam_files(bam_file_path)
  })
  output$bam_files <- renderUI({
    checkboxGroupInput("select_bam_files", label = h4("Select the Relevant 
                                                      Bam FIles"),
                       choices = found_bam_files(), 
                       selected = found_bam_files()[1])   
  })
  processed_gff <- reactive({
    withProgress(message = 'Processing the gff file.',
                 detail = 'This only happens once', value = 0,{
      modify_gff_inplace(paste(gff_file_path, input$gff_select,sep=""))
                 }
    )
  })
  output$gff_check <- renderPrint({
    cat("The gff file may take a moment to be processed, the plot will appear when finished")
  })
  
  gffInput <- reactive({
  filter_gff_for_rows(processed_gff(), input$select_genes)
  })
  output$gff_rows<- renderDataTable({
    gffInput () 
  })

    select_group_fun <- reactive({
      count <- 1
      lapply(find_bam_files(bam_file_path)
        [1:length(input$select_bam_files)], 
        function(i) {          
            selectInput(paste0('snumber', i),              
            h5(paste0('Select a group for ', i)),
            choices = 1:length(input$select_bam_files))
        }
      )
    })
 output$select_group <- renderUI({
   select_group_fun()
  })    
 group_list <- reactive({
   res <- lapply(find_bam_files(bam_file_path)
                 [1:length(input$select_bam_files)], 
                 function(i) { 
                   input[[paste0('snumber', i)]]
                 })
 })
 
  
  poly_a_counts<- reactive({
 
    initial_table <- get_a_counts (bam_file_path, gffInput(),input$select_bam_files,
                                   input$select_genes,group_list())
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
    make_plot(poly_a_counts(), input$xslider,input$select_genes, input$legend, input$merge)
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

  
})