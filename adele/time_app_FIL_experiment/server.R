library(shiny)
#library(nesoni)
library(varistran)
library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)
load("proprocessed_counts.RDS")

#options(bitmapType='cairo')
shinyServer(function(input, output) {
  
  
  dat_count <- reactive({
    withProgress(message = "Processing data", value = 0, {
      
      #Filter  all the columns by the user input
      dat_count <- count_mat_preprocessing[apply(count_mat_preprocessing, 1, function(row) {any(row >= input$num_filter)}),]
      
      #Use varistran to normalise
      dat_count <- varistran::vst(dat_count, cpm = T)
      
      # Pull out the names of the columns
      all_names <- colnames(dat_count)
      
      # Create a vector of the YPD columns and timepoints 7 and 8
      excluded_names <- c("B1_YPD_R1", "B1_YPD_R2", "B2_YPD_R1", "B2_YPD_R2", "R1_FIL_7", "R2_FIL_7", "R1_FIL_8", "R2_FIL_8",
                          "R1_INHIB_7", "R2_INHIB_7", "R1_INHIB_8", "R2_INHIB_8")
      
      # Subset down to the columns of interest
      kept_names <- all_names[! all_names %in% excluded_names]
      
      # Remove the YPD columns
      dat_count <- dat_count[, kept_names]
      
      return(dat_count)
    }) 
  })
  
  #Generate a selectInput UI element that contains the names of genes
  output$gene_select <- renderUI({
    
    # Pull out the rownames from dat_count
    gene_list <- rownames(dat_count())
    
    # Create a named vector
    names(gene_list) <- gene_list
    
    # Create a selectizeInput for the gene names and pre-load a gene of interest ("ARG8")
    selectizeInput("gene_selected", label = h4("Find a gene of interest"),
                   choices = gene_list,
                   selected = gene_list["ARG8/orf19.3770"], multiple = FALSE, options = list(maxOptions = 5))
    
  })
  
  # Get the summary statistics (maximum and minimum expression across all genes at each time point) for the second view mode
  summ_stats <- reactive({
    
    # Get the counts data frame
    data <- dat_count()
    
    # Create a vector of the time points used in the experiment (e.g 2HR=2.00, 2HR15MIN=2.25, 2HR30MIN=2.50, 2HR45MIN=2.75, 3HR=3.00, 4HR=4.00)
    Time <- as.numeric(c(2.00, 2.15, 2.30, 2.45, 3.00, 4.00))
    
    # Create two empty vectors and an empty data frame
    maximums <- c()
    minimums <- c()
    df <- data.frame()
    
    # Loop through the timepoints, grab out all the columns at a given timepoint and create a data-frame for that time point 
    # (e.g time point 1 (2 hours) includes R1_FIL_1, R2_FIL_1, R1_INHIB_1, R2_INHIB_1. These are spaced out by 6 columns in the counts data frame). 
    # Populate the empty maximum/minimum vectors by taking the max/min of the dataframe at the given timepoint
    for(i in 1:length(Time)){
      col1 <- data[, i, drop = F]
      col2 <- data[, i + 6, drop = F]
      col3 <- data[, i + 12, drop = F]
      col4 <- data[, i + 18, drop = F]
      df <- data.frame(col1, col2, col3, col4)
      maximums[i] <- max(df)
      minimums[i] <- min(df)
    }
    
    # Bind the maximum and minimum values for the different time points into a dataframe
    summ_stats <- data.frame(
      Time,
      maximums,
      minimums
    )
    return(summ_stats)
  })
  
  # Convert the counts data frame from wide format to long. This works better for plotting
  melted_count <- reactive({
    
    # Pull in the counts data-frame
    melted_count <- dat_count()
    
    # Grab the column names from the counts data frame
    nam <- colnames(melted_count)
    
    # Create a vector of the time points used in the experiment (e.g 2HR=2.00, 2HR15MIN=2.25, 2HR30MIN=2.50, 2HR45MIN=2.75, 3HR=3.00, 4HR=4.00)
    time_points <- as.numeric(c(2.00, 2.25, 2.50, 2.75, 3.00, 4.00))
    
    # Repeat the vector four times (this is because the column names are ordered sequentially)
    all_time <- rep(time_points, 4)
    
    # Create an empty vector
    new_names <- c()
    
    # Loop through the column names. Pull out the replicate and condition type, then paste together with the new time point
    # i.e this step is just to change the columns names (from something like: R1_FIL_1, R1_FIL_2, R1_FIL_3 to R1_FIL_2, R1_FIL_2.25, R1_FIL_2.5)
    for (x in 1:length(nam)) {
      d <- unlist(strsplit(nam[x], "_"))[[1]]
      e <- unlist(strsplit(nam[x], "_"))[[2]]
      new_names[x] <- paste0(d, "_", e, "_", all_time[x])
    }
    
    # Give the counts matrix the new column names
    colnames(melted_count) <- new_names
    
    # Melt the counts matrix
    melted_count <- melt(as.matrix(melted_count))
    
    # Fix the column names
    colnames(melted_count)[1] <- "Gene"
    colnames(melted_count)[3] <- "Expression"
    
    # Separate the second column into three columns, the replicate number (R1/R2), the condition (INHIB/FIL) and which time point
    data <- melted_count %>%
      separate(Var2, into = c("Replicate", "Condition", "Time"), sep = "_")
    
    return(data)
  })
  
  
  ########################################################
  ########################################################
  #And now for outputs
  
  #Generate the download button
  output$download_eps <- downloadHandler(
    filename = function(){
      paste0("mDivi_plot_", Sys.Date(), ".eps")
    },
    content = function(file) {
      ggsave(file)
    })
  
  
  #Generate the line plot
  output$time_plot <- renderPlot({
    withProgress(message = "Plotting data...", value = 0, {
      
      # Get the melted count data frame and the summary statistics
      df <- melted_count()
      stats <- summ_stats()
      
      # Create the base layer of the plot

      m <- ggplot(data = stats, aes(x = Time)) + scale_x_discrete(name = "Time (HR)", limits = c(2.00, 2.25, 2.50, 2.75, 3.00, 4.00),
                                                                  labels = c(2.00, 2.15, 2.30, 2.45, 3.00, 4.00))
      
      # If the first view mode is selected, just line plot the gene of interest. If the second mode is selected, add a ribbon for the summary statistics
      if(input$plot_mode == 1) {
        m <- m + geom_line(data = subset(df, Gene == input$gene_selected), aes(x = as.numeric(Time),
                                                                               y = Expression, group = 
                                                                                 interaction(Replicate, Condition),
                                                                               colour = Condition, linetype = Replicate), size = 1.5) +
          ylab("Expression (log2 CPM)") + theme_bw() 

      } else if (input$plot_mode == 2) {
        m <- m + geom_ribbon(aes(ymin = minimums, ymax = maximums), colour = "#cccccc", alpha = 0.2)
        m <- m + geom_line(data = subset(df, Gene == input$gene_selected), aes(x = as.numeric(Time),
                                                                               y = Expression, 
                                                                               group = interaction(Replicate, Condition),
                                                                               colour = Condition, 
                                                                               linetype = Replicate), size = 1.5) +
          ylab("Expression (log2 CPM)") + theme_bw()

       } 
      return(m)
    })
  })

  #Give plot info
  output$info_txt <- renderText({

    paste("time=", input$plot_click$x, "   expression =", input$plot_click$y)

  })
  

})

