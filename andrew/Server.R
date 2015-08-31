# a change
source("major_calculations.R")
library("Rsamtools")
library("jsonlite")

shinyServer(function(input, output, session) {
  output$select_file_path <- renderUI({
    selectInput("file_path", label = h4("Select a dataset"), 
                choices = list.dirs(full.names=F, recursive =F), 
                selected =  list.dirs(full.names=F, recursive =F)[1])   
  })
  
  found_gff_files <- reactive({
    find_gff_files(paste0(input$file_path))
  })
  output$gff_files <- renderUI({
    selectInput("gff_select", label = h4("GFF File Selection"), 
                choices = found_gff_files(), 
                selected =  found_gff_files()[1])   
  })
  found_bam_files <- reactive({
    find_bam_files(paste0(input$file_path, "/"))
  })
  output$bam_files <- renderUI({
    
    if (class(found_bam_files())=='data.frame'){
      checkboxGroupInput("select_bam_files", label = h4("Select the Relevant 
                                                        Bam FIles"),
                         choices = found_bam_files()[[1]], 
                         selected = found_bam_files()[[1]][1]) 
      #Goes into the data frame and gets the file paths corresponding 
      #to the selected BAM files      
    }
    else{
      checkboxGroupInput("select_bam_files", label = h4("Select the Relevant 
                                                        Bam FIles"),
                         choices = found_bam_files(), 
                         selected = found_bam_files()[1])  
    }
    
  })
  processed_gff <- reactive({
    withProgress(message = 'Processing the gff file, this may take a few seconds.',
                 detail = 'This only happens once.', value = 0,{   
                   if (substring(input$gff_select,1,1)=="/"){
                     modify_gff_inplace(paste(input$gff_select))
                   }
                   else{
                     modify_gff_inplace(paste("./", input$file_path,"/", input$gff_select,sep=""))
                   }
                   
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
    lapply(input$select_bam_files, 
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
    res <- lapply(input$select_bam_files, 
                  function(i) { 
                    input[[paste0('snumber', i)]]
                  })
  })
  
  
  poly_a_counts<- reactive({
    if (class(found_bam_files())=='data.frame'){
      bam_files <- found_bam_files()$bam [found_bam_files()$name ==input$select_bam_files]
    }
    else{
      bam_files <- input$select_bam_files
    }
    withProgress(message = 'Processing the bam files.',
                 detail = 'This will take longer if there are a lot of reads.', 
                 value = 0,{  
                   
                   initial_table <- get_a_counts (input$file_path, gffInput(),bam_files,
                                                  group_list(),found_bam_files())
                   initial_table[is.na(initial_table)]<-0
                   
                   if (input$all_reads ==F){  
                     initial_table <- initial_table[initial_table$number_of_as > 0,]
                   }
                   subsetted_by_sliders <- initial_table[initial_table$number_of_ad_bases >=            
                                                           input$ad_slider&
                                                           initial_table$width >=
                                                           input$al_length[1]&
                                                           initial_table$width <=
                                                           input$al_length[2],]                        
                 })
  })
  output$means_frame <- renderDataTable({
    make_means_and_meds_frame(poly_a_counts())
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
    names_string (split_frame, input$merge, input$all_reads) 
  })
  plot_calcs <- reactive({
    make_plot(poly_a_counts(), input$xslider,input$select_genes, input$legend, 
              input$merge, input$alt_plot, input$order_alt, input$alt_cumu_dis, input$poly_a_pileup)
  })    
  output$scp_plot<- renderPlot({  
    plot_calcs()
  })
  #Workaround for a shiny bug thatdoesn't handle reactive plots well. 
  plot_calcs2 <- function(){
    make_plot(poly_a_counts(), input$xslider,input$select_genes, input$legend, 
              input$merge, input$alt_plot, input$order_alt, input$alt_cumu_dis, input$poly_a_pileup)
  }
  selected_plot_points <- reactive({
    input$plot_brush
  })
  output$seq_plot <- renderDataTable({
    min_points <-selected_plot_points()$xmin
    max_points <-selected_plot_points()$xmax
    print(str(poly_a_counts()))
    if (input$alt_plot ==F){
      sequence <- poly_a_counts()$sequence[poly_a_counts()$number_of_as > min_points & poly_a_counts()$number_of_as < max_points]  
      name <- poly_a_counts()$gene_or_peak_name[poly_a_counts()$number_of_as > min_points & poly_a_counts()$number_of_as < max_points] 
      sample <- poly_a_counts()$sample[poly_a_counts()$number_of_as > min_points & poly_a_counts()$number_of_as < max_points]
      
    }
    else{
      sequence <- poly_a_counts()$sequence[poly_a_counts()$width > min_points & poly_a_counts()$width < max_points]  
      name <- poly_a_counts()$gene_or_peak_name[poly_a_counts()$width > min_points & poly_a_counts()$width < max_points] 
      sample <- poly_a_counts()$sample[poly_a_counts()$width > min_points & poly_a_counts()$width < max_points]
    }
    data.frame(sample, name, sequence)
  })
  
  output$frame_type <-renderText({
    if (input$alt_plot ==F){
      return("Reads That Fall Within The Selected Poly (A)-Tail Lengths")
    }
    else{
      return("Reads That Fall Within The Selected Read Lengths")      
    }
  }) 

  output$downloadPlot <- downloadHandler(
    filename = function(){
      paste(trim(input$select_genes), '.eps', sep='')
    },
    content = function(file){
      setEPS(width = 10)
      postscript(file)
      plot_calcs2()     
      dev.off()       
    })  
})