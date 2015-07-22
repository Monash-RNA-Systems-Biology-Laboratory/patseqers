#Server.R
#Authors - Michael See
#        - Paul Harrison
library(shiny)
library(limma)
library(edgeR)
library(nesoni)
library(DT)

source("help.R")

shinyServer(function(input,output) {
    
    output$selDataSet <- renderUI({        
        sets.df <- read.csv("setname.csv", header=T)
        sets.df <- data.frame(lapply(sets.df, as.character), stringsAsFactors=F)
        stvec <- as.character(sets.df$setname)
        selectInput("dataset", "Choose dataset",sets.df$setname)
        
    })
    
    output$selDataSet1 <- renderUI({
        sets.df <- read.csv("setname.csv", header=T)
        sets.df <- data.frame(lapply(sets.df, as.character), stringsAsFactors=F)
        stvec <- as.character(sets.df$setname)
        selectInput("dataset1", "Choose dataset",sets.df$setname)
    })
    
    output$selCol <- renderUI({
        name.l <- colnames(x1())
        selectizeInput("col2Disp", "Select", choices = name.l, multiple = T, selected=NULL)
    })
    
    x1 <- reactive({
        if(input$tab == 0) {
            tmp <- input$dataset
        } else {
            tmp <- input$dataset1
        }
        x1.df <- data.frame(read.csv(file.path(
            paste0(tmp),
            paste0(tmp,"-count.csv")), 
            stringsAsFactors=F, 
            header = TRUE))
        
        return(x1.df)
    })
    
    y1 <- reactive({
        if(input$tab == 0) {
            tmp <- input$dataset
        } else {
            tmp <- input$dataset1
        }
        
        y1.df <- data.frame(read.csv(file.path(
            paste0(tmp),
            paste0(tmp,"-info.csv")), 
            stringsAsFactors=F, 
            header = TRUE))
        
        return(y1.df)
    })
    
    minN <- reactive({
        return(input$nMin)
    })
    spanN <- reactive({
        if(input$tab == 0){
            return(input$minSpan)
        } else {
            #return(4.0)
            return(input$minSpan1)
        }
    })
    priorC <-  reactive({
        return(input$prc)
    })
    
    flt <-reactive({
        hld <- x1()
        rownames(hld) <- hld$Name 
        hld <- hld[,-1]
        flt.df <-data.frame(hld[apply(hld,1,function(row) {
            any(row >= input$nMin)
        }),])
        
        return(flt.df)
    })
#     output$plotui <- renderUI({
#         plotOutput("plot1", width="1400", height="850",
#                    
#                    brush = brushOpts(
#                        id = "plot_brush",                       
#                        direction = 'y',#input$brush_dir,
#                        resetOnNew = T
#                    )
#         )
#     })
    mrg <- reactive({
        x2 <- flt()
        x2 <- DGEList(x2)
        x2 <- calcNormFactors(x2)
        x2 <- data.frame(cpm(x2, log=T, prior.count=priorC()))
        
        x2$Name <- rownames(x2)
        y2 <- merge(x2, y1())
        t1 <- data.frame(y2$gene, y2$product)
        colnames(t1) <- c("Gene", "Product")
        return(t1)
    })
    
    
    pMake <- reactive({       
        x2 <- flt()
        x2 <- DGEList(x2)
        x2 <- calcNormFactors(x2)
        x2 <- data.frame(cpm(x2, log=T, prior.count=priorC()))
        if(length(input$col2Disp==0)){
            x2 <- x2[,-which(names(x2) %in% input$col2Disp)]
        }
        x2$Name <- rownames(x2)
        y2 <- merge(x2, y1())
        
        xy.ls <- list()
        xy.ls$E <- as.matrix(subset(x2, select=-c(Name)))
        xy.ls$genes <- data.frame(y2$gene, y2$product)
        colnames(xy.ls$genes) <- c("gene", "product")
        rownames(xy.ls$genes) <- x2$Name
        xy.elist <- new("EList", xy.ls)   
        if(input$tab==1){
            xy.ls$genes <- xy.ls$genes[input$datab_rows_all,]
            xy.ls$E <- xy.ls$E[input$datab_rows_all,]
            xy.elist <- new("EList", xy.ls)
        }
        if(length(input$datab_rows_selected) != 0){
            xy.ls$genes <- xy.ls$genes[input$datab_rows_selected,]
            xy.ls$E <- xy.ls$E[input$datab_rows_selected,]
            xy.elist <- new("EList", xy.ls)
        }
        a <- t1("prefix", xy.elist, min.span=spanN())
        return(a)
    })
    
    
    output$plot1 <- renderPlot({
        a <- pMake()
        return(a)
    })
    
    output$datab = DT::renderDataTable(mrg(), server = F, rownames=F,
                                       options = list(searchHighlight = TRUE
                                                      ))

    output$tout = renderPrint({
        cat('\n\nFiltered Rows\n\n')
        #cat(input$datab_rows_all, sep = ', ')
        cat(length(input$datab_rows_all))
        cat('\n\nSelected rows:\n\n')
        #cat(input$datab_rows_selected, sep = ', ')
        cat(length(input$datab_rows_selected))
    })


    
}    
)
