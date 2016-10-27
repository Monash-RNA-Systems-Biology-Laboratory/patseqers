library(edgeR)
library(biomaRt)
library(ggplot2)
library(nesoni)
source("helper.R")

all_folders <- list.dirs(full.names=F, recursive =F)
selected_directory <- all_folders[!grepl("www", all_folders)]

file_call <- find_files(selected_directory)
#peak_file_call <- find_peaks(selected_directory)

selected.mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "celegans_gene_ensembl", host = "www.ensembl.org")
go.attribute <- "refseq_mrna"

shinyServer(function(input, output) {
  
  
  #Select which organism database from ensembl biomart will use
  # selected.mart <- reactive({
  #   mart <- "celegans_gene_ensembl"
  #   mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = mart, host = "www.ensembl.org")
  #   return(mart)
  # })
  # 
  # go.attribute <- reactive({
  #   
  #   attr <- "refseq_mrna"
  #   return(attr)
  # })
  
  
  #Create a reactive object by pulling out the counts dataframe and filtering by the input$numfilter
  count.df <- reactive({
    count.df <- file_call[[1]]
    
    xsel <- input$x.selection
    ysel <- input$y.selection
    sel <- c(xsel, ysel)
    count.df <- count.df[, sel, drop = FALSE]
    if(input$type_filter == 1) {
      count.df <- count.df[apply(count.df, 1, function(row) {any(row >= input$numfilter)}), ]
    } else if(input$type_filter == 2) {
      count.df <- count.df[apply(count.df, 1, function(row) {all(row >= input$numfilter)}), ]
    }
    count.df <- varistran::vst(count.df, cpm = T)
    return(count.df)
  })
  
  #Create reactive objects by pulling out the tail length and tail counts dataframe
  
  tail.count.df <- reactive({
    tail.count.df <- file_call[[3]]
    xsel <- input$x.selection
    ysel <- input$y.selection
    sel <- c(xsel, ysel)
    tail.count.df <- tail.count.df[, sel, drop = FALSE]
    return(tail.count.df)
  })
  
  #Filter the tail length file by the tail.counts.df
  tail.df <- reactive({
    tail.df <- file_call[[2]]
    tail.count <- tail.count.df()
    
    xsel <- input$x.selection
    ysel <- input$y.selection
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
    xsel <- input$x.selection
    
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
    ysel <- input$y.selection
    
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
    #files <- file.call()
    #files <- file_call
    info.df <- file_call[[4]]
    
    info.df$significant <- rownames(info.df) %in% rownames(base.df)
    info.df <- subset(info.df, significant)
    info.df$significant <- NULL
    return(info.df)
  })
  
  base.df <- reactive({
    withProgress(message = "Building dataframe with selected inputs", 
                 value = 0, {
                   base.df <- initial.df()
                   info.df <- info.df()
                   base.df$gene <- info.df$gene
                   return(base.df)
                 })
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
  output$gene.search.ui <- renderUI({
    withProgress(message = "Finding genes...", value = 0, { 
      selectizeInput("gene.search.sel", label = "Input genes to highlight", choices = gene.names(), 
                     selected = NULL, multiple = TRUE, options = list(maxOptions = 10))
    })
  })
  
  
  
  #Create a dataframe filtered by the selected genes. This will take input from both the text input
  #and the selectiveInput widget
  selected.genes <- reactive({
    
    if(!is.null(input$gene.search.sel) | !is.null(input$gene.paste.sel)){
      base.df <- base.df()
      rownames(base.df) <- tolower(rownames(base.df))
      base.df$gene <- tolower(base.df$gene)
      
      sel.genes <- tolower(input$gene.search.sel)
      pas.genes <- tolower(input$gene.paste.sel)
      
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
    
    goterm.word <- input$goterm.select
    goslim.word <- input$goslim.select
    org <- selected.mart#()
    attr <- go.attribute#()
    
    
    if(input$type_goterm == 1  & !goterm.word == "") {
      query <- fetch_goterm(attr, goterm.word, org)
      colnames(query) <- "genes"
      
      
      if(nrow(query) != 0) {
        base.df <- base.df()
        base.df <- base.df[unique(grep(paste(query$genes, collapse = "|"), rownames(base.df),
                                       value = TRUE, ignore.case = TRUE)), ]
        return(base.df)
      } 
    } else if (input$type_goterm == 2 & !goslim.word == "") {
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
                 value = 0, {
                   df <- base.df()
                   
                   geneplot <- ggplot(df, aes(x= x.mean, y = y.mean)) + geom_point() + theme_bw() + 
                     xlab(paste("x.mean", input$datatypex)) + #theme(aspect.ratio =1) +
                     ylab(paste("y.mean", input$datatypey)) #+  xlim(vmin, vmax) + ylim(vmin, vmax)
                   
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
                   
                   
                   if(!input$goterm.select == "" | input$type_goterm == 2  & !input$goslim.select == ""){ 
                     goterm.df <- goterm.genes()
                     if(is.null(goterm.df) == FALSE){
                       geneplot <- geneplot + geom_point(data = goterm.df, aes(x = x.mean, y = y.mean), shape = 21,
                                                         colour = "#3366ff", fill = "#3399ff", size = 2.5)
                     }
                   } else {
                     geneplot
                   }
                   
                   if(!is.null(input$gene.search.sel) | !is.null(input$gene.paste.sel)){
                     sel.genes.df <- selected.genes()
                     geneplot <- geneplot + geom_point(data = sel.genes.df, aes(x = x.mean, y=y.mean), colour = "#ff9900",
                                                       fill = "#ffff00", shape = 23, size = 2.5)
                   } else{
                     geneplot
                   }
                   
                   
                   
                   
                   geneplot
                 })
  })
  
  output$gplot <- renderPlot({
    
    main_plot()
  })
  
  #Reports on the number of genes in total being plotted
  output$total.txt <- renderText({
    total.genes <- nrow(base.df())
    
    return(paste(total.genes, "total genes in the graph after filtering. Below, the top 10 closest genes to
the clicked mouse position, arranged in order of nearest distance to furthest. To find a gene more accurately, use the Search Genes tab."))
  })
  
  #Returns the closest rows to the mouse hover
  output$info_plot <- renderPrint({
    df <- base.df()
    rownames(df) <- NULL
    df <- df[,c(3, 1, 2)]
    nearPoints(df, input$plot_click, xvar = "x.mean", yvar = "y.mean", threshold = 10, 
               maxpoints = 10)
  })
  
  #Reports on the number of genes highlighted by the gene search widgets
  output$sel.txt <- renderText({
    sel.genes <- nrow(selected.genes())
    
    if(!is.null(input$gene.search.sel) | !is.null(input$gene.paste.sel)){
      
      return(paste(sel.genes, "genes have been selected and highlighted in the current dataset."))
    } 
    
  })
  
  #Gives the output for the selected genes in a table.
  output$sel.table <- renderTable({
    
    if(!is.null(input$gene.search.sel) | !is.null(input$gene.paste.sel)){
      df <- selected.genes()
      rownames(df) <- NULL
      df$selected <- NULL
      return(df)
    } 
  })
  
  
  
  #Reports the number of genes matching the GO Term
  output$goterm.txt <- renderText({
    
    
    if(!input$goterm.select == "" | input$type_goterm == 2  & !input$goslim.select == ""){
      df <- goterm.genes()
      rownames(df) <- NULL
      sel.goterm <- nrow(goterm.genes())
      if(is.null(goterm.genes()) == FALSE){
        return(paste(sel.goterm, "genes have been annotated with the selected GO term in the current 
                         dataset."))
      } else {
        return(paste("0 genes have been annotated with the selected GO term in the current 
                         dataset."))
      }
      
    } else {
      return(paste("0 genes have been annotated with the selected GO term in the current 
                         dataset."))
    }
  })
  
  #Gives the output for the GO Term selected genes in a table.
  output$goterm.table <- renderTable({
    if(!input$goterm.select == "" | input$type_goterm == 2 & !input$goslim.select == "" ){
      df <- goterm.genes()
      rownames(df) <- NULL
      
     # if(is.null(df == FALSE)){
        return(df)
    #  }
    }
  })
  
  output$deps <- downloadHandler(
    filename = function(){
      paste0(input$file_name, ".eps")
    },
    content = function(file) {
      ggsave(file)
    })
  
  output$dpdf <- downloadHandler(
    filename = function(){
      paste0(input$file_name, ".pdf")
    },
    content = function(file) {
      ggsave(file)
    })
  
  output$dsel <- downloadHandler(
    filename = function(){
      paste0(input$file_name, "selected.csv")
    },
    content = function(file) {
      write.csv(selected.genes(), file)
    })
  
  
  output$dgot <- downloadHandler(
    filename = function(){
      paste0(input$file_name, "goterm.csv")
    },
    content = function(file) {
      write.csv(goterm.genes(), file)
    })
  
  
  
})

