library(edgeR)
library(biomaRt)
library(ggplot2)
library(shinyBS)
#library(nesoni)
#library(dplyr)
#library(tidyr)
source("helper.R")

# Load the RDS object
load("www/data/preload_files_called.rds")

#Load the mart object
load("www/data/preload_mart.rds")

shinyServer(function(input, output, session) {
  
  
  xselection <- reactive({
    
    
    
    xsel <- switch(input$input_x,
                   "N2" = c("N2.rep1", "N2.rep2", "N2.rep3"),
                   "gld2" = c("Gld2.rep1", "Gld2.rep2", "Gld2.rep3"),
                   "gld2_pcb19_RNAi" = c("Gld2.pcb19.rep1", "Gld2.pcb19.rep2"),
                   "gld2_ccf1_RNAi" = c("Gld2.CCF1.rep1", "Gld2.CCF1.rep2"),
                   "1.2_cell" = c("X1.2.cell.egg.rep1", "X1.2.cell.egg.rep2")
    )
    return(xsel)
  })
  
  yselection <- reactive({
    
    ysel <- switch(input$input_y,
                   "N2" = c("N2.rep1", "N2.rep2", "N2.rep3"),
                   "gld2" = c("Gld2.rep1", "Gld2.rep2", "Gld2.rep3"),
                   "gld2_pcb19_RNAi" = c("Gld2.pcb19.rep1", "Gld2.pcb19.rep2"),
                   "gld2_ccf1_RNAi" = c("Gld2.CCF1.rep1", "Gld2.CCF1.rep2"),
                   "1.2_cell" = c("X1.2.cell.egg.rep1", "X1.2.cell.egg.rep2")
    )
    return(ysel)
  })
  
  #Create a reactive object by pulling out the counts dataframe and filtering by the input$numfilter
  count.df <- reactive({
    count.df <- file_call[[1]]
    
    xsel <- xselection()
    ysel <- yselection()
    sel <- c(xsel, ysel)
    
    count.df <- count.df[, sel, drop = FALSE]
    if(input$type_filter == 1) {
      count.df <- count.df[apply(count.df, 1, function(row) {any(row >= input$numfilter)}), ]
    } else if(input$type_filter == 2) {
      count.df <- count.df[apply(count.df, 1, function(row) {all(row >= input$numfilter)}), ]
    }
    #count.df <- varistran::vst(count.df, cpm = T)
    return(count.df)
  })
  
  #Create reactive objects by pulling out the tail length and tail counts dataframe
  
  tail.count.df <- reactive({
    tail.count.df <- file_call[[3]]
    xsel <- xselection()
    ysel <- yselection()   
    sel <- c(xsel, ysel)
    
    tail.count.df <- tail.count.df[, sel, drop = FALSE]
    return(tail.count.df)
  })
  
  #Filter the tail length file by the tail.counts.df
  tail.df <- reactive({
    tail.df <- file_call[[2]]
    tail.count <- tail.count.df()
    
    
    xsel <- xselection()
    ysel <- yselection()
    sel <- c(xsel, ysel)
    
    tail.df <- tail.df[, sel, drop = FALSE]
    if(input$type_filter == 1) {
      tail.df <- tail.df[apply(tail.count, 1, function(row) {any(row >= input$numfilter)}), ]
    } else if(input$type_filter == 2) {
      tail.df <- tail.df[apply(tail.count, 1, function(row) {all(row >= input$numfilter)}), ]
    }
    return(tail.df)
  })
  
  #Generate a dataframe containing the columns picked for X and y axis from either the counts or tail length 
  select.xaxis.rep <- reactive({
    xsel <- xselection()
    
    if(input$datatypex == "gene expression") {
      count.df <- count.df()
      x.df <- count.df[, xsel, drop = FALSE]
      x.df <- tmmnorm(x.df)
    } else if(input$datatypex == "tail length") {
      tail.df <- tail.df()
      x.df <- tail.df[, xsel, drop = FALSE]
    } 
    return(as.data.frame(x.df))
  })
  
  
  select.yaxis.rep <- reactive({
    ysel <- yselection()
    
    if(input$datatypey == "gene expression") {
      count.df <- count.df()
      y.df <- count.df[, ysel, drop = FALSE]
      y.df <- tmmnorm(y.df)
    } else if(input$datatypey == "tail length"){
      tail.df <- tail.df()          
      y.df <- tail.df[, ysel, drop = FALSE]
    } 
    return(as.data.frame(y.df))
  })
  
  #Find the means for the selected replicates on the xaxis and yaxis
  xaxis.mean <- reactive({
    df <- select.xaxis.rep()
    x.mean <- rowMeans(df)
    x.mean <- as.data.frame(x.mean)
    return(x.mean)
  })
  
  yaxis.mean <- reactive({
    df <- select.yaxis.rep()
    y.mean <- rowMeans(df)
    y.mean <- as.data.frame(y.mean)
    return(y.mean)
  })
  
  #Generate a dataframe holding the means of the samples selected on the x and y axis
  #If the xaxis is different from the yaxis (ie if one is gene expression data and the other is tail 
  #length), then they will have been filtered differently and the resulting dataframes will have a different 
  #number of rows. The genes in common must be found and then the dataframes for each are subsetted
  #to the genes in common before being added together in one dataframe
  
  initial.df <- reactive({
    xaxis <- xaxis.mean()
    yaxis <- yaxis.mean()
    
    xaxis$significant <- rownames(xaxis) %in% rownames(yaxis)
    yaxis$significant <- rownames(yaxis) %in% rownames(xaxis)
    
    xaxis <- subset(xaxis, significant)
    xaxis$significant <- NULL
    yaxis <- subset(yaxis, significant)
    yaxis$significant <- NULL
    
    df <- data.frame(xaxis, yaxis)
    df <- df[!(is.na(df$x.mean) | is.na(df$y.mean)), ]
    
    return(df)
  })
  
  
  
  
  #Create a reactive object by pulling out the info dataframe and subset it by the base dataframe to ensure
  #that it matches the results of the current filtering. It can't be filtered by the filter number from UI
  #as the info.df has no count data. This way of filtering ensures that it matches regardless of tail 
  #count or gene count filtering
  
  info.df <- reactive({
    base.df <- initial.df()
    info.df <- file_call[[4]]
    
    info.df$significant <- rownames(info.df) %in% rownames(base.df)
    info.df <- subset(info.df, significant)
    info.df$significant <- NULL
    return(info.df)
  })
  
  
  #From the info.df, both the locus names (rownames of info.df) and the common/sequence names
  #can be retrieved from the genes column in info.df (column 2 usually). This can be used to create
  #a list of genes for a selectizeinput element in the ui that'll allow genes to be selected.
  
  gene.names <- reactive({
    info.df <- info.df()
    gname <- as.character(info.df[,2])
    names <- c(gname)
    names <- as.list(names)
    
    return(names)
  })
  
  #Create a ui element for the list of genes currently present in the plot
  output$gene_search_ui <- renderUI({
     
      selectizeInput("gene_search_ui", label = "Input genes to highlight", choices = gene.names(), 
                     selected = NULL, multiple = TRUE, options = list(maxOptions = 10), width = "80%")
  })
  
  addTooltip(session, id = "gene_search_ui", title = 
               "Help: Only 10 gene names will be displayed at a time but you can type into the box to find the gene you are interested it. If you can't find the gene you are looking for, it may not be in the dataset or the filter may need to be adjusted",
             placement = "top", trigger = "hover")
  
  ### Create a base dataframe 
  base.df <- reactive({
    withProgress(message = "Building dataframe with selected inputs", 
                 value = 0, style = "old", {
                   
                   base.df <- initial.df()
                   info.df <- info.df()
                   base.df$gene <- info.df$gene
                   return(base.df)
                 })
  })
  
  
  #Create a dataframe filtered by the selected genes. This will take input from both the text input
  #and the selectiveInput widget
  selected.genes <- reactive({
    
    if(!is.null(input$gene_search_ui) | !is.null(input$gene_paste_sel)){
      base.df <- base.df()
      rownames(base.df) <- tolower(rownames(base.df))
      base.df$gene <- tolower(base.df$gene)
      
      sel.genes <- tolower(input$gene_search_ui)
      pas.genes <- tolower(input$gene_paste_sel)
      
      sgenes <- as.data.frame(sel.genes)
      
      pgenes <- strsplit(pas.genes, " ")[[1]]
      pgenes <- as.data.frame(pgenes)
      
      base.df$selected <- rownames(base.df) %in% sgenes$sel.genes | base.df$gene %in% sgenes$sel.genes |
        rownames(base.df) %in% pgenes$pgenes | base.df$gene %in% pgenes$pgenes
      base.df <- subset(base.df, selected)
      return(base.df)
    }
    
  })
  
  
  #Search by GoTerms. fetch_goterm is defined in helper.R and it runs a biomaRt query dependent on three 
  #arguements. 
  
  goterm.genes <- reactive({
    
    goterm.word <- input$goterm_select
    goslim.word <- input$goslim_select
    org <- selected.mart
    attr <- "refseq_mrna"
    
    
    if(input$type_goterm == 2  & !goterm.word == "") {
      query <- fetch_goterm(attr, goterm.word, org)
      colnames(query) <- "genes"
      
      
      if(nrow(query) != 0) {
        base.df <- base.df()
        base.df <- base.df[unique(grep(paste(query$genes, collapse = "|"), rownames(base.df),
                                       value = TRUE, ignore.case = TRUE)), ]
        return(base.df)
      } 
    } else if (input$type_goterm == 1 & !goslim.word == "") {
      query <- fetch_goterm(attr, goslim.word, org)
      colnames(query) <- "genes"
      
      
      if(nrow(query) != 0) {
        base.df <- base.df()
        base.df <- base.df[unique(grep(paste(query$genes, collapse = "|"), rownames(base.df),
                                       value = TRUE, ignore.case = TRUE)), ]
        return(base.df)
      } 
    }
  })
  
  
  
  
  ####Below here are outputs#########################
  
  #Renders the plot
  main_plot <- reactive({
    withProgress(message = "Generating graph...", 
                 value = 0, style = "old", {
                   df <- base.df()
                   
                   geneplot <- ggplot(df, aes(x= x.mean, y = y.mean)) + geom_point() +
                     xlab(paste(input$input_x, input$datatypex)) + ylab(paste(input$input_y, input$datatypey)) +
                     theme_bw() + theme(aspect.ratio=1)
                   # + 
                   #theme(aspect.ratio=1)
                   #+  xlim(vmin, vmax) + ylim(vmin, vmax)
                   #theme(aspect.ratio =1) +
                   #               if(input$plot_axes == 2) {
                   if(input$datatypex == input$datatypey) {
                     
                     xmean.max <- max(df$x.mean)
                     ymean.max <- max(df$y.mean)
                     xmean.min <- min(df$x.mean)
                     ymean.min <- min(df$y.mean)
                     vmax <- max(c(xmean.max, ymean.max))
                     vmin <- min(c(xmean.min, ymean.min))
                     
                     geneplot <- geneplot + xlim(vmin, vmax) + ylim(vmin, vmax)
                     print(geneplot)
                   } else {
                     geneplot
                   }
                   
                   
                   if(!input$goterm_select == "" | input$type_goterm == 1  & !input$goslim_select == ""){ 
                     goterm.df <- goterm.genes()
                     if(is.null(goterm.df) == FALSE){
                       geneplot <- geneplot + geom_point(data = goterm.df, 
                                                         aes(x = x.mean, y = y.mean, 
                                                             fill = "GO Term Annotated genes"), 
                                                         shape = 21, size = 2.5)
                     }
                   } else {
                     geneplot
                   }
                   
                   if(!is.null(input$gene_search_ui) | !is.null(input$gene_paste_sel)){
                     sel.genes.df <- selected.genes()
                     if(nrow(sel.genes.df) >= 1)  {
                       geneplot <- geneplot + geom_point(data = sel.genes.df, aes(x = x.mean, y=y.mean, fill = "Individual genes"), 
                                                         shape = 23, size = 2.5)
                     }
                   } else{
                     geneplot
                   }
                   
                   
                   geneplot <- geneplot + scale_fill_discrete(name = "Selected data")
                   geneplot
                 })
  })
  
  
  # Render the plot
  output$gplot <- renderPlot({
    
    main_plot()
    
  })
  
  # Return a table summarising the number of genes found
  output$report_stats <- renderTable({
    total.genes <- nrow(base.df())
    
    
    if(!is.null(input$gene_search_ui) | !is.null(input$gene_paste_sel)){
      sel.genes <- nrow(selected.genes())
      
    } else {
      sel.genes <- as.integer(0)
    }
    
    if(!input$goterm_select == "" | input$type_goterm == 1  & !input$goslim_select == ""){
      if(is.null(goterm.genes()) == FALSE){
        sel.goterm <- nrow(goterm.genes())
        
      } else {
        sel.goterm <- as.integer(0)
      }
      
    } else {
      sel.goterm <- as.integer(0)
    }
    
    report_df <- data.frame(
      total.genes,
      sel.genes,
      sel.goterm
    )
    colnames(report_df) <- c("Total genes", "Individual genes selected", "Genes annotated with selected GO Term/Slim")
    
    return(report_df)
  })
  
  
  #Returns the closest rows to the mouse hover
  output$info_plot <- renderPrint({
    df <- base.df()
    rownames(df) <- NULL
    df <- df[,c(3, 1, 2)]
    nearPoints(df, input$plot_click, xvar = "x.mean", yvar = "y.mean", threshold = 10, 
               maxpoints = 5)
  })
  
  
  # Gives the output for the individually selected genes in a table.
  output$sel.table <- renderTable({
    
    if(!is.null(input$gene_search_ui) | !is.null(input$gene_paste_sel)){
      df <- selected.genes()
      rownames(df) <- NULL
      df$selected <- NULL
      return(df)
    } 
  })
  
  # Returns a table of the genes that have been annotated with the GO Term or Go Slim of interest
  output$annotated_selected <- renderTable({
    
    df <- goterm.genes()
    return(df)
  })
  
  # Save figure as a EPS file
  output$deps <- downloadHandler(
    filename = function(){
      paste0(input$file_name, ".eps")
    },
    content = function(file) {
      ggsave(file)
    })
  
  # Save figure as a PDF file
  output$dpdf <- downloadHandler(
    filename = function(){
      paste0(input$file_name, ".pdf")
    },
    content = function(file) {
      ggsave(file)
    })
  # Save selected genes into a csv
  output$dsel <- downloadHandler(
    filename = function(){
      paste0(input$file_name, "selected.csv")
    },
    content = function(file) {
      write.csv(selected.genes(), file)
    })
  
  
  # Save GO Term/Slim annotated genes into a csv 
  output$dgot <- downloadHandler(
    filename = function(){
      paste0(input$file_name, "goterm.csv")
    },
    content = function(file) {
      write.csv(goterm.genes(), file)
    })
  
})

