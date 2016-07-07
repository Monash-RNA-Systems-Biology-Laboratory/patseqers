library(shiny)
#library(nesoni)
library(varistran)
library(dplyr)
library(reshape2)
library(tidyr)
#library(limma)
#library(edgeR)
#library(Matrix)
#library(parallel)
library(ggplot2)
source("helper.R")
#options(bitmapType='cairo')
shinyServer(function(input, output) {
  
  #Find the counts.csv file within the tail-tools pipeline output
  #Counts.csv will be a grouped table
  file_call <- reactive({
    files <- find_file(input$select_experiment)
    return(files)
  })
  
  dat_count <- reactive({
    withProgress(message = "Retrieving files", value = 0, {
      #Retrieve the reactive grouped table objec
      df_count <- file_call()
      
      #Pull out the count matrix and the annotations for the count matrix
      dat_count <- df_count$Count
      annot_df <- df_count$Annotation
      
      annot_ensem <- rownames(annot_df)
      annot_gene <- annot_df$gene
      
      gene_names <- c()
      for(i in 1:length(annot_ensem)) {
        gene_names[i] <- paste0(annot_gene[i], "_", annot_ensem[i])
      }
      
      rownames(dat_count) <- gene_names
      
      
      #Filter  all the columns by the user input
      #dat_count <- dat_count[apply(dat_count, 1, function(row) {any(row >= input$num_filter)}),]
      
      #Use varistran to normalise
      dat_count <- varistran::vst(dat_count, cpm = T)
      #print(attr(dat_count, "lib.size"))
      return(dat_count)
    }) 
  })
  
  #Generate a selectInput UI element that contains the names of genes
  output$gene_select <- renderUI({
    
    gene_list <- rownames(dat_count())
    selectizeInput("gene_selected", label = h4("Select gene"),
                   choices = gene_list,
                   selected = gene_list, multiple = FALSE, options = list(maxOptions = 5))
    
  })
  
  
  main_names <- reactive({
    df_names <- colnames(dat_count())
    excluded_names <- c("YPD_R1", "YPD_R2")
    df_names <- df_names[! df_names %in% excluded_names]
    return(df_names)
  })
  
  summ_stats <- reactive({
    data <- dat_count()[, main_names()]
    Time <- as.numeric(c(2.00, 2.15, 2.30, 2.45, 3.00, 4.00, 5.00, 6.00))
    
    maximums <- c()
    minimums <- c()
    df <- data.frame()
    
    for(i in 1:length(Time)){
      df <- data[, i, drop = F]
      maximums[i] <- max(df)
      minimums[i] <- min(df)
      
    }
    
    summ_stats <- data.frame(
      Time,
      maximums,
      minimums
    )
    return(summ_stats)
  })
  
  
  melted_count <- reactive({
    melted_count <- dat_count()[, main_names()]
    nam <- colnames(melted_count)
    time_points <- as.numeric(c(2.00, 2.25, 2.50, 2.75, 3.00, 4.00, 5.00, 6.00))
    all_time <- rep(time_points, 2)
    new_names <- c()
    
    for (x in 1:length(nam)) {
      d <- unlist(strsplit(nam[x], "_"))[[1]]
      new_names[x] <- paste0(d, "_", all_time[x])
    }
    colnames(melted_count) <- new_names
    
    melted_count <- melt(as.matrix(melted_count))
    colnames(melted_count)[1] <- "Gene"
    colnames(melted_count)[3] <- "Expression"
    data <- melted_count %>%
      separate(Var2, into = c("Condition", "Time"), sep = "_")
    return(data)
  })
  
  melted_ypd <- reactive({
    dat <- dat_count()[, c(1:2)]
    mean_dat <- as.data.frame(rowMeans(dat), drop = F)
    colnames(mean_dat) <- "YPD"
    
    melted_dat <- melt(as.matrix(mean_dat))
    colnames(melted_dat) <- c("Gene", "Condition", "Expression")
    return(melted_dat)
  })
  
  test_ypd <- reactive({
    dat <- dat_count()[, c(1:2)]
    
    melted_dat <- melt(as.matrix(dat))
    colnames(melted_dat) <- c("Gene", "Condition", "Expression")
    return(melted_dat)
  })
  
  
  ########################################################
  ########################################################
  #And now for outputs
  
  #Generate the download button
  output$download_eps <- downloadHandler(
    filename = function(){
      paste0(input$sel.gene, Sys.Date(), ".eps")
    },
    content = function(file) {
      ggsave(file)
    })
  
  
  #Generate the line plot
  output$time_plot <- renderPlot({
    withProgress(message = "Loading data...", value = 0, {
      df <- melted_count()
      dat_ypd <- melted_ypd()
      stats <- summ_stats()
      test_dat <- test_ypd()
      
      m <- ggplot(data = stats, aes(x = Time))
      
      if(input$plot_mode == 1) {
        m <- m + geom_line(data = subset(df, Gene == input$gene_selected), aes(x = as.numeric(Time),
                                                                               y = Expression, group = Condition,
                                                                               colour = Condition), size = 1.5) + 
          scale_x_discrete(name = "Time", limits = c(2.00, 2.25, 2.50, 2.75, 3.00, 4.00, 5.00, 6.00),
                           labels = c(2.00, 2.15, 2.30, 2.45, 3.00, 4.00, 5.00, 6.00)) + theme_bw()
        
      } else if (input$plot_mode == 2) {
        m <- m + geom_ribbon(aes(ymin = minimums, ymax = maximums), colour = "#cccccc", alpha = 0.2)
        m <- m + geom_line(data = subset(df, Gene == input$gene_selected), aes(x = as.numeric(Time),
                                                                               y = Expression, group = Condition,
                                                                               colour = Condition), size = 1.5) +
          scale_x_discrete(name = "Time", limits = c(2.00, 2.25, 2.50, 2.75, 3.00, 4.00, 5.00, 6.00),
                           labels = c(2.00, 2.15, 2.30, 2.45, 3.00, 4.00, 5.00, 6.00)) + theme_bw()
        
      } else {
        m <-  m + geom_line(data = subset(df, Gene == input$gene_selected), aes(x = as.numeric(Time),
                                                                                y = Expression, group = Condition,
                                                                                colour = Condition), size = 1.5) +
          geom_hline(data = subset(test_dat, Gene ==input$gene_selected), aes(yintercept = Expression, group = Condition, colour = Condition), size = 1.5) +
          scale_x_discrete(name = "Time", limits = c(2.00, 2.25, 2.50, 2.75, 3.00, 4.00, 5.00, 6.00),
                           labels = c(2.00, 2.15, 2.30, 2.45, 3.00, 4.00, 5.00, 6.00)) + theme_bw()
      }
      
      
      return(m)
    })
  })
  
  #Give plot info
  output$info_txt <- renderText({
    
    paste("time=", input$plot_click$x, "   expression =", input$plot_click$y)
    
  })
  
  # Check what data is actually being plotted
  output$table_plotted <- renderTable({
    print(subset(melted_count(), Gene ==input$gene_selected))
  })
  
  #Check the library size
  output$library_size_table <- renderTable({
    dt <- data.frame(true_library_size = as.integer(attr(dat_count(), "true.lib.size")),
      adjusted_library_size = as.integer(attr(dat_count(), "lib.size")))
    print(dt)
  })
  
})

