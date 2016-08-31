library(edgeR)
library(biomaRt)
library(ggplot2)
library(nesoni)
source("helper.R")

selected_directory <- list.dirs(full.names=F, recursive =F)
file_call <- find_files(selected_directory)
peak_file_call <- find_peaks(selected_directory)

shinyServer(function(input, output) {
  
  #Create a reactive object containing the selected directory
  # selected.directory <- reactive({
  #     return(list.dirs(full.names=F, recursive =F))
  # })  
  
  #Select which organism database from ensembl biomart will use
  selected.mart <- reactive({
    mart <- "celegans_gene_ensembl"
    mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = mart, host = "www.ensembl.org")
    return(mart)
  })
  
  go.attribute <- reactive({
    
    attr <- "refseq_mrna"
    return(attr)
  })
  
  #Use the find_files function to find the following files in the selected directory: "genewise-count.csv", 
  #"genewise-tail.csv and "genewise-info.csv". These will be loaded in the format of a list containing
  #three dataframes. The first item in the list is the counts dataframe. 2nd item is the tail dataframe
  #and the third is the info dataframe.
  # file.call <- reactive({
  #     withProgress(message = "Finding files", value = 0, {
  #         files <- find_files(selected.directory())
  #         return(files)
  #     })
  # })
  
  #Create a reactive object containing column names 
  
  axis.names <- reactive({
    
    # files <- file.call()
    files <- file_call
    count.df <- files[[1]]
    axis.names <- colnames(count.df)
    return(as.list(axis.names))
  })
  
  #Use the axis.names list to render selectInput in UI for x and y axis
  output$x.sel.ui <- renderUI({
    
    selectInput("x.selection", label = h4("Select samples for the x-axis"),
                choices = axis.names(),
                selected =  axis.names()[1], multiple = TRUE)  
    
  })
  
  output$y.sel.ui <- renderUI({
    
    selectInput("y.selection", label = h4("Select samples for the y-axis"),
                choices = axis.names(),
                selected =  axis.names()[2], multiple = TRUE)  
    
  })
  
  
  
  #Create a reactive object by pulling out the counts dataframe and filtering by the input$numfilter
  count.df <- reactive({
    # files <- file.call()
    
    files <- file_call
    count.df <- files[[1]]
    
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
    # files <- file.call()
    files <- file_call
    tail.count.df <- files[[3]]
    xsel <- input$x.selection
    ysel <- input$y.selection
    sel <- c(xsel, ysel)
    tail.count.df <- tail.count.df[, sel, drop = FALSE]
    return(tail.count.df)
  })
  
  #Filter the tail length file by the tail.counts.df
  tail.df <- reactive({
    # files <- file.call()
    files <- file_call
    tail.df <- files[[2]]
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
  
  #Create a reactive object containing the peaks    
  peak.df <- reactive({
    # df <- find_peaks(selected.directory())
    # tail.df <- df$Tail_count        
    tail.df <- peak_file_call$Tail_count
    xsel <- input$x.selection
    ysel <- input$y.selection
    sel <- c(xsel, ysel)
    
    pro.df <- peak_file_call$Proportion
    
    pro.df <- pro.df[,grep(paste(sel, collapse = "|"), colnames(pro.df), ignore.case = T)]
    if(input$type_filter == 1) {
      pro.df <- pro.df[apply(tail.df, 1, function(row) {any(row >= input$numfilter)}), ]
      
    } else {
      pro.df <- pro.df[apply(tail.df, 1, function(row) {all(row >= input$numfilter)}), ]
    }        
    peak1 <- pro.df[, grep("peak1", colnames(pro.df))]
    peak2 <- pro.df[, grep("peak2", colnames(pro.df))]
    final <- rownames(pro.df)
    for(i in 1:ncol(peak1)){
      final <- data.frame(final, peak1[,i] - peak2[,i])
    }
    colnames(final) <- c("gene", colnames(count.df()))
    
    names <- as.character(final$gene)
    for(i in 1:length(names)){
      names[i] <- unlist(strsplit(names[i], "-peak"))[[1]]
    }
    
    rownames(final) <- names    
    final$gene <- NULL
    return(final)
  })
  
  #Unnecessary reactive objects, can just use input$x.selection directly in the select.xaxis.rep(licates)
  #However, may have used the xselect objects elsewhere, need to ensure this before deleting from code 
  #completely
  
  #Generate a dataframe containing the columns picked for X and y axis from either the counts or tail length 
  select.xaxis.rep <- reactive({
    xsel <- input$x.selection
    
    if(input$datatypex == "gene expression") {
      count.df <- count.df()
      x.df <- count.df[, xsel, drop = FALSE]
      # x.df <- tmmnorm(x.df)
    } else if(input$datatypex == "tail length") {
      tail.df <- tail.df()
      x.df <- tail.df[, xsel, drop = FALSE]
    } else {
      peak.df <- peak.df()
      x.df <- peak.df[, xsel, drop = FALSE]
    }
    return(as.data.frame(x.df))
  })
  
  
  select.yaxis.rep <- reactive({
    ysel <- input$y.selection
    
    if(input$datatypey == "gene expression") {
      count.df <- count.df()
      y.df <- count.df[, ysel, drop = FALSE]
      
      #y.df <- tmmnorm(y.df)
    } else if(input$datatypey == "tail length"){
      tail.df <- tail.df()          
      y.df <- tail.df[, ysel, drop = FALSE]
    } else {
      peak.df <- peak.df()
      y.df <- peak.df[, ysel, drop = FALSE]
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
    files <- file_call
    info.df <- files[[4]]
    
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
    #         rname <- row.names(info.df)
    gname <- as.character(info.df[,2])
    #         names <- c(gname, rname)
    names <- c(gname)
    names <- as.list(names)
    
    return(names)
  })
  
  #Create a ui element for the list of genes currently present in the plot
  output$gene.search.ui <- renderUI({
    selectizeInput("gene.search.sel", label = "Input genes to highlight", choices = gene.names(), 
                   selected = NULL, multiple = TRUE, options = list(maxOptions = 5))
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
  
  #Generate a dataframe containing genes matching a key word in the product description. The product 
  #description is found in info.df 
  key.genes <- reactive({
    
    if(!input$key.select == ""){
      base.df <- base.df()
      key.word <- input$key.select
      info.df <- info.df()
      
      base.df$key_search <- info.df$product
      base.df <- base.df[grep(key.word, base.df$key_search, ignore.case = TRUE), ]
      return(base.df)
    } 
    
    
  })
  
  #Search by GoTerms. fetch_goterm is defined in helper.R and it runs a biomaRt query dependent on three 
  #arguements. 
  
  goterm.genes <- reactive({
    
    goterm.word <- input$goterm.select
    org <- selected.mart()
    attr <- go.attribute()
    
    if(!goterm.word == ""){
      query <- fetch_goterm(attr, goterm.word, org)
      colnames(query) <- "genes"
      
      
      if(nrow(query) != 0) {
        base.df <- base.df()
        base.df <- base.df[unique(grep(paste(query$genes, collapse = "|"), rownames(base.df),
                                       value = TRUE, ignore.case = TRUE)), ]
        #          base.df <- base.df[unique(grep(paste(query$genes, collapse = "|"), rownames(base.df),
        #                                                value = TRUE, ignore.case = TRUE)), ]
        #          base.df <- base.df[grep(paste(query$genes, collapse = "|"), rownames(base.df),
        #                                            value = TRUE, ignore.case = TRUE), ]
        # print(head(query))
        # print(str(query))
        # print(length(query))
        # print(nrow(query))
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
                   
                   if(input$plot_axes == 2) {
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
                   
                   
                   if(!is.null(input$gene.search.sel) | !is.null(input$gene.paste.sel)){
                     sel.genes.df <- selected.genes()
                     geneplot <- geneplot + geom_point(data = sel.genes.df, aes(x = x.mean, y=y.mean), colour = "#ff00ff",
                                                       fill = "#b366ff", shape = 23)
                   } else{
                     geneplot
                   }
                   # 
                   #                          if(!input$key.select == "") {
                   #                              key.genes.df <- key.genes()
                   #                              geneplot <- geneplot + geom_point(data = key.genes.df, aes(x = x.mean, y=y.mean), 
                   #                                                                colour = "#FFF000")
                   #                          } else {
                   #                              geneplot
                   #                          }
                   
                   
                   if(!input$goterm.select == ""){
                     goterm.df <- goterm.genes()
                     geneplot <- geneplot + geom_point(data = goterm.df, aes(x = x.mean, y = y.mean), 
                                                       colour = "#3399ff")
                   } else {
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
    
    return(paste(total.genes, "total genes in the graph after filtering. Below, the top 5 closest genes to
the clicked mouse position, arranged in order of nearest distance to furthest. To find a gene more accurately, use the Data Select tab."))
  })
  
  #Returns the closest rows to the mouse hover
  output$info_plot <- renderPrint({
    df <- base.df()
    rownames(df) <- NULL
    df <- df[,c(3, 1, 2)]
    nearPoints(df, input$plot_click, xvar = "x.mean", yvar = "y.mean", threshold = 10, 
               maxpoints = 5)
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
      
      return(selected.genes())
    } 
    
  })
  
  #Reports the number of genes found matching the key word
  output$key.txt <- renderText({
    sel.key <- nrow(key.genes())
    
    if(!input$key.select == "") {
      return(paste(sel.key, "genes match the key word searched for in the current dataset."))
    }
  })
  
  #Gives the output for the key word selected genes in a table.
  output$key.table <- renderTable({
    if(!input$key.select == "") {
      
      return(key.genes())
    } 
    #         return(head(initial.df()))
  })
  
  #Reports the number of genes matching the GO Term
  output$goterm.txt <- renderText({
    
    
    if(!input$goterm.select == ""){
      sel.goterm <- nrow(goterm.genes())
      return(paste(sel.goterm, "genes have been annotated with the selected GO term in the current 
                         dataset."))
    } 
  })
  
  #Gives the output for the GO Term selected genes in a table.
  output$goterm.table <- renderTable({
    if(!input$goterm.select == ""){
      
      return(goterm.genes())
    } 
    #         return(head(base.df()))
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
  
  output$dkey <- downloadHandler(
    filename = function(){
      paste0(input$file_name, "keyterm.csv")
    },
    content = function(file) {
      write.csv(key.genes(), file)
    })
  
  output$dgot <- downloadHandler(
    filename = function(){
      paste0(input$file_name, "goterm.csv")
    },
    content = function(file) {
      write.csv(goterm.genes(), file)
    })
  
  
  
})

