library(shiny)
library(ggplot2)
library(varistran)
library(reshape2)
library(nesoni)
library(dplyr)
library(tidyr)
source("helper.R")

shinyServer(function(input, output){
    
    #Find the selected directory file path
    selected.directory <- reactive({
        return(input$select.experiment)
    })
    
    #Use the directory file path to find the counts.csv file - csv contains a grouped table
    file.call <- reactive({
        files <- find_file(selected.directory())
        return(files)
    })
    
    #Pull the count data out of the grouped table, filter the values and then normalise with varistran
    dat_count <- reactive({
        
        #Retrieve the grouped table reactive object
        df_count <- file.call()
        
        #Pull out the count matrix and the annotations for the count matrix
        dat_count <- df_count$Count
        
        annot <- df_count$Annotation
        annot_ensem <- rownames(annot)
        annot_gene <- annot$gene
        
        gene_names <- c()
        for(i in 1:length(annot_ensem)) {
            gene_names[i] <- paste0(annot_gene[i], "_", annot_ensem[i])
        }
        
        rownames(dat_count) <- gene_names        
        
        
        #Filter  all the columns by the user input
        dat_count <- dat_count[apply(dat_count, 1, function(row) {all(row >= input$num.filter)}),]
        
        #Use varistran to normalise
        dat_count <- varistran::vst(dat_count)
        
        return(dat_count)    
    })
    
    gene_list <- reactive({
        dat_count <- dat_count()
        gene_list <- rownames(dat_count)
        return(as.list(gene_list))
    })
    
    output$gene.select <- renderUI({
        
        selectizeInput("sel.gene", label = h4("Select gene"),
                    choices = gene_list(),
                    selected =  gene_list()[1], multiple = FALSE, options = list(maxOptions = 10))  
        
    })
    
    new_names <- reactive({
        df_names <- colnames(dat_count())
        
        samples_removed <- c()
        sample_sec <- c()
        time_sec <- c()
        replicate_sec <- c()
        
        for(i in 1:length(df_names)) {
            samples_removed[i] <- unlist(strsplit(df_names[i], split = "Sample_", fixed = TRUE))[[2]]
        }
        
        for(i in 1:length(df_names)) {
            sample_sec[i] <- substr(samples_removed[i], 4, nchar(samples_removed[i])-1)
            time_sec[i] <- substr(samples_removed[i], nchar(samples_removed[i]), nchar(samples_removed[i]))
            replicate_sec[i] <- substr(samples_removed[i], 1, 2)
        }
        
        new_names <- c()
        for(i in 1:length(df_names)) {
            new_names[i] <- paste0(sample_sec[i], "_", time_sec[i], "_", replicate_sec[i])            
        }
        return(new_names)
    })
    
    replicate_names <- reactive({
        if(input$org.sel == 1) {
            replicate_names <- c("Ca_1_E", "Ca_6_E", "Ca_9_E", "CaMac_1_E", "CaMac_6_E", "CaMac_9_E")
        } else {
            replicate_names <- c("Mac_1_E", "Mac_6_E", "Mac_9_E", "CaMac_1_E", "CaMac_6_E", "CaMac_9_E")
        }
        return(replicate_names)
    })
    
    mean_count <- reactive({
        dat_count <- dat_count()
        new_names <- new_names()
        replicate_names <- replicate_names()
        colnames(dat_count) <- new_names
        
        mean_dat <- data.frame()
        
        for(i in 1:length(replicate_names)) {
            if(i == 1) {
                mean_dat <- rowMeans(dat_count[, grep(paste0("^", replicate_names[i], "."), colnames(dat_count))])
            } else {
                mean_dat <- cbind(mean_dat, rowMeans(dat_count[, grep(paste0("^", replicate_names[i], "."),
                                                                      colnames(dat_count))]))
            }
        }
        mean_names <- c()
        for(i in 1:length(replicate_names)) {
            mean_names[i] <- substr(replicate_names[i], 1, nchar(replicate_names[i]) - 2)
        }
        colnames(mean_dat) <- mean_names
        
        return(mean_dat)
    })
    
    summ_stats <- reactive({
        data <- mean_count()
        timepoints <- as.factor(c(1, 6, 9))
        
        maximuns <- c()
        minimuns <- c()
        df <- data.frame()
        
        for(i in 1:length(timepoints)){
            df <- data[, grep(paste0("_", timepointst[i]), colnames(data))]
            maximums[i] <- max(df[,1:2])
            minimums[i] <- min(df[,1:2])
            
        }
        
        summ_stats <- data.frame(
            timepoints,
            maximums,
            minimums
            )
        
    })
    
    melted_count <- reactive({
        mean_count <- mean_count()
        melted_count <- melt(mean_count)
        
        colnames(melted_count)[1] <- "Gene"
        colnames(melted_count)[3] <- "Expression"
        data <- melted_count %>%
            separate(Var2, into = c("Condition", "Time"), sep = "\\_")
        
        return(data)
    })
    
    #Outputs
    output$time.plot <- renderPlot({
        df <- melted_count()
        ggplot(data = subset(df, Gene == input$sel.gene), aes(as.factor(Time), Expression, 
                                                              group = Condition, colour = Condition)) + 
            geom_line() + theme_bw() + theme(aspect.ratio =1) + xlab("Timepoints (hr)") + 
            ylab("Gene expression (log2 normalised counts)")
#         ggplot(data = df, aes(as.factor(Time)))
    })
    
    output$info.txt <- renderText({
        
        paste("time=", input$plot_click$x, "   expression =", input$plot_click$y)
        
    })
})