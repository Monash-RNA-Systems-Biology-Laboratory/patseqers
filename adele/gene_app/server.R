library(edgeR)
library(biomaRt)
library(shinyBS)
library(ggplot2)
library(nesoni)
source("helper.R")

shinyServer(function(input, output) {
    # setwd("/data/home/apattison/ShinyApps/andrew/")
    output$experi_sel <- renderUI({
        selectInput("select.experiment", label = "Select experiment", 
                    choices = list.dirs(full.names=F, recursive =F), 
                    selected =  list.dirs(full.names=F, recursive =F)[1])   
    })
    
    #Create a reactive object containing the selected directory
    selected.directory <- reactive({
        return(input$select.experiment)
    })  
    
    #Select which organism database from ensembl biomart will use
    selected.mart <- reactive({
        mart <- input$org.type
        mart <- useMart("ensembl", dataset = mart)
        return(mart)
    })
    
    go.attribute <- reactive({
        org <- input$org.type
        if(input$org.type == "hsapiens_gene_ensembl"){
            attr <- "ucsc"
        } else if(input$org.type == "scerevisiae_gene_ensembl"){
            attr <- "ensembl_gene_id"
        } else if(input$org.type == "celegans_gene_ensembl") {
            attr <- "refseq_mrna"
        }
        return(attr)
    })
    
    #Use the find_files function to find the following files in the selected directory: "genewise-count.csv", 
    #"genewise-tail.csv and "genewise-ifno.csv". These will be loaded in the format of a list containing
    #three dataframes. The first item in the list is the counts dataframe. 2nd item is the tail dataframe
    #and the third is the info dataframe.
    file.call <- reactive({
        files <- find_files(selected.directory())
        return(files)
    })
    
    #Create a reactive object by pulling out the counts dataframe and filtering by the input$numfilter
    count.df <- reactive({
        files <- file.call()
        count.df <- files[[1]]
        
        count.df <- count.df[apply(count.df, 1, function(row) {any(row >= input$numfilter)}), ]
        return(count.df)
    })
    
    #Create reactive objects by pulling out the tail length and tail counts dataframe
    
    tail.count.df <- reactive({
        files <- file.call()
        tail.count.df <- files[[3]]
        return(tail.count.df)
    })
    
    #Filter the tail length file by the tail.counts.df
    tail.df <- reactive({
        files <- file.call()
        tail.df <- files[[2]]
        tail.count <- tail.count.df()
        tail.df <- tail.df[apply(tail.count, 1, function(row) {any(row >= input$numfilter)}), ]   
        return(tail.df)
    })
    
    
    #Create a reactive object containing column names from count.df
    
    axis.names <- reactive({
        axis.names <- colnames(count.df())
        return(as.list(axis.names))
    })
    
    #Use the axis.names list to render selectInput in UI for x and y axis
    output$x.sel.ui <- renderUI({
        
        selectInput("x.selection", label = h4("Select X-axis"),
                    choices = axis.names(),
                    selected =  axis.names()[1], multiple = TRUE)  
        
    })
    
    output$y.sel.ui <- renderUI({
        
        selectInput("y.selection", label = h4("Select Y-axis"),
                    choices = axis.names(),
                    selected =  axis.names()[2], multiple = TRUE)  
        
    })
    
    peak.df <- reactive({
        df <- find_peaks(selected.directory())
        tail.df <- df$Tail_count        
        
        pro.df <- df$Proportion
        pro.df <- pro.df[apply(tail.df, 1, function(row) {any(row >= input$numfilter)}), ]
        
        peak1 <- pro.df[, grep("peak1", colnames(pro.df))]
        peak2 <- pro.df[, grep("peak2", colnames(pro.df))]
        final <- rownames(pro.df)
        for(i in 1:ncol(peak1)){
            final <- data.frame(final, peak1[,1] - peak2[,i])
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
        count.df <- count.df()
        tail.df <- tail.df()
        peak.df <- peak.df()
        xsel <- input$x.selection
        #         xsel <- xselect()
        if(input$datatypex == "gene expression") {
            x.df <- count.df[, xsel, drop = FALSE]
            x.df <- tmmnorm(x.df)
        } else if(input$datatypex == "tail length") {
            x.df <- tail.df[, xsel, drop = FALSE]
        } else {
            x.df <- peak.df[, xsel, drop = FALSE]
        }
        return(as.data.frame(x.df))
    })
    
    
    select.yaxis.rep <- reactive({
        count.df <- count.df()
        tail.df <- tail.df()
        peak.df <- peak.df()
        ysel <- input$y.selection
        #         ysel <- yselect()
        
        if(input$datatypey == "gene expression") {
            y.df <- count.df[, ysel, drop = FALSE]
            y.df <- tmmnorm(y.df)
        } else if(input$datatypey == "tail length"){
            y.df <- tail.df[, ysel, drop = FALSE]
        } else {
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
    
    base.df <- reactive({
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
        base.df <- base.df()
        files <- file.call()
        info.df <- files[[4]]
        
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
            sel.genes <- tolower(input$gene.search.sel)
            pas.genes <- tolower(input$gene.paste.sel)
            info.df <- info.df()
            base.df$genes <- tolower(info.df$gene)
            
            sgenes <- as.data.frame(sel.genes)
            
            pgenes <- strsplit(pas.genes, " ")[[1]]
            pgenes <- as.data.frame(pgenes)
            
            base.df$selected <- rownames(base.df) %in% sgenes$sel.genes | base.df$genes %in% sgenes$sel.genes |
                rownames(base.df) %in% pgenes$pgenes | base.df$genes %in% pgenes$pgenes
            base.df <- subset(base.df, selected)
            return(base.df)
        }
        
    })
    
    #Generate a dataframe containing genes matching a key word in the product description. The product 
    #description is found in info.df 
    key.genes <- reactive({
        base.df <- base.df()
        if(!input$key.select == ""){
            
            key.word <- input$key.select
            info.df <- info.df()
            base.df$key_search <- info.df$product
            base.df <- base.df[grep(key.word, base.df$key_search, ignore.case = TRUE), ]
        } 
        
        return(base.df)
    })
    
    #Search by GoTerms. fetch_goterm is defined in helper.R and it runs a biomaRt query dependent on three 
    #arguements. 
    
    goterm.genes <- reactive({
        base.df <- base.df()
        goterm.word <- input$goterm.select
        org <- selected.mart()
        attr <- go.attribute()
        
        if(!goterm.word == ""){
            query <- fetch_goterm(attr, org, goterm.word)
            colnames(query) <- "genes"
            base.df <- base.df[unique(grep(paste(query$genes, collapse = "|"), rownames(base.df),
                                           value = TRUE, ignore.case = TRUE)), ]
            #             base.df$goterm <- rownames(base.df) %in% query$genes
            #             base.df <- subset(base.df, goterm)
            #             base.df$goterm <- NULL
        }
        return(base.df)
    })
    
    ####Below here are outputs#########################
    
   
    
    #     m.df <- reactive({
    #         df <- base.df()
    #         x <- input$plot_hover$x
    #         y <- input$plot_hover$y
    #         
    #         
    #         m.df <- df[df$x.mean >= (x - 0.5) & df$x.mean <= (x + 0.5) &
    #                        df$y.mean >= (y - 0.5) & df$y.mean <= (y + 0.5), ]
    #         print(paste(nrow(m.df), "X min =", x - 0.5, "X max =", x + 0.5,
    #                     "Y min =", y - 0.5, "Y max =", y + 0.5), str(m.df))
    #         return(m.df)
    #     })
    #     
    #     output$hover.txt <- renderTable({
    #         df <- m.df()
    #         print(paste(nrow(m.df), "X min =", x - 0.5, "X max =", x + 0.5,
    #                     "Y min =", y - 0.5, "Y max =", y + 0.5), str(m.df))
    #         return(head(df))
    #     })
    
    
    
    #Renders the plot
    main_plot <- reactive({
        df <- base.df()
        sel.genes.df <- selected.genes()
        goterm.df <- goterm.genes()
        key.genes.df <- key.genes()
        xmean.max <- max(df$x.mean)
        ymean.max <- max(df$y.mean)
        xmean.min <- min(df$x.mean)
        ymean.min <- min(df$y.mean)
        vmax <- max(c(xmean.max, ymean.max))
        vmin <- min(c(xmean.min, ymean.min))
        
        geneplot <- ggplot(df, aes(x= x.mean, y = y.mean)) + geom_point() + theme_bw() + 
            xlab(paste("x.mean", input$datatypex)) + #theme(aspect.ratio =1) +
            ylab(paste("y.mean", input$datatypey)) + xlim(vmin, vmax) + ylim(vmin, vmax)
        
        geneplot
        
        
        if(!is.null(input$gene.search.sel) | !is.null(input$gene.paste.sel)){
            
            geneplot <- geneplot + geom_point(data = sel.genes.df, aes(x = x.mean, y=y.mean), 
                                              colour = "#FF00FF")
        } else{
            geneplot
        }
        
        if(!input$key.select == "") {
            geneplot <- geneplot + geom_point(data = key.genes.df, aes(x = x.mean, y=y.mean), 
                                              colour = "#FFF000")
        } else {
            geneplot
        }
        
        
        if(!input$goterm.select == ""){
            
            geneplot <- geneplot + geom_point(data = goterm.df, aes(x = x.mean, y = y.mean), 
                                              colour = "#0000FF")
        } else {
            geneplot
        }
        
        
        geneplot
    })
    
    output$gplot <- renderPlot({
        
        main_plot()
    })
    
    #Reports on the number of genes in total being plotted
    output$total.txt <- renderText({
        total.genes <- nrow(base.df())
        
        return(paste(total.genes, "total genes in the graph after filtering. Below, the top 5 closest genes to
the current mouse position, arranged in order of nearest distance to furthest."))
    })
    
    #Returns the closest rows to the mouse hover
    output$info_plot <- renderPrint({
        
        nearPoints(base.df(), input$plot_hover, xvar = "x.mean", yvar = "y.mean", threshold = 10, 
                   maxpoints = 5)
    })
    
    #Reports on the number of genes highlighted by the gene search widgets
    output$sel.txt <- renderText({
        sel.genes <- nrow(selected.genes())
        
        if(!is.null(input$gene.search.sel) | !is.null(input$gene.paste.sel)){
        
            return(paste(sel.genes, "genes have been selected and highlighted in the current dataset."))
        } 
        
    })
    
    #Reports the number of genes found matching the key word
    output$key.txt <- renderText({
        sel.key <- nrow(key.genes())
        
        if(!input$key.select == "") {
            return(paste(sel.key, "genes match the key word searched for in the current dataset."))
        }
    })
    
    
    #Reports the number of genes matching the GO Term
    output$goterm.txt <- renderText({
        sel.goterm <- nrow(goterm.genes())
        
        if(!input$goterm.select == ""){
            
            return(paste(sel.goterm, "genes have been annotated with the selected GO term in the current 
                         dataset."))
        } 
    })
    
    output$deps <- downloadHandler(
        filename = function(){
            paste0(input$file_name, ".eps")
        },
        content = function(file) {
            setEPS(width = 10)
            postscript(file)
            print(main_plot())
            dev.off()
        })
    
    output$dpdf <- downloadHandler(
        filename = function(){
            paste0(input$file_name, ".pdf")
        },
        content = function(file) {
            ggsave(file)
        })
    
    
})

