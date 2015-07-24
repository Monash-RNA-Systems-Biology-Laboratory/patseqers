#Server.R
#Authors - Michael See    mmsee2@student.monash.edu.au
#        - Paul Harrison

# A slightly long winded side-note on brush interactivity in this app:
# Brushes (At least the way I have enabled them) for one reason or another does not have the resolution to 
# select all the rows in the heatmap. 
# To get around this, one may consider grabbing the index of the first and last rows the brush returns 
# and match that to the dataframe and pull out all the rows in between to match and pull out the rows it needs.
# Keeping in mind however, that the brush also fails to accurately select the last row brushed as well, so 
# realistically, interactivity with the dendogram is might be more accurate if one worked out how to enable 
# interactivity in that region of the plot. This is all of course something to consider while keeping in mind the 
# plot is not interactive on the right region, and only allows brushes to be drawon over the rowmeans column.



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
        selectizeInput("col2Disp", "Select a column to not display", choices = name.l, multiple = T, selected=NULL)
    })
    output$selCol1 <- renderUI({
        name.l <- colnames(x1())
        selectizeInput("col2Disp1", "Select a column to not display", choices = name.l, multiple = T, selected=NULL)
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
    output$plotui <- renderUI({
        input$gopt
        plotOutput("plot1", width=isolate(input$pwidth), height=isolate(input$pheight)
# Commented out brush functions as they are a little broken and does not select all the rows being brushed

#                    brush = brushOpts(
#                        id = "plot_brush",                       
#                        direction = 'y',#input$brush_dir,
#                        resetOnNew = F
#                   )
        
        )
    })

    # Called when the datatable needs to be rendered
    # Contains the information that goes into that table
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
    
    #Function which plots
    #Great assistance from Paul in debugging some side effects of reactive expressions
    pMake <- function()({
        
        x2 <- flt()
        
        #Two if statements de-select columns for plotting
        #Would only need one if shiny would realise one is not being rendered by a conditional panel
        if(length(input$col2Disp==0)){
            x2 <- x2[,-which(names(x2) %in% input$col2Disp)]
        }
        if(length(input$col2Disp1==0)){
            x2 <- x2[,-which(names(x2) %in% input$col2Disp1)]
        }
        
        x2 <- DGEList(x2)
        x2 <- calcNormFactors(x2)
        x2 <- data.frame(cpm(x2, log=T, prior.count=priorC()))
        
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
        
        #t1 now returns a dataframe which is set to xy, but this is not currently used
        xy <<- t1("prefix", xy.elist, min.span=spanN()) 
            
    })
    
    pht <- reactive({
        return(input$pheight)
    })
    pwt <- reactive({
        return(input$pwidth)
    })
    
    output$plot1 <- renderPlot({

        a <- pMake()
        return(a)
    })
    #Renders the datatable
    #The pMake() function takes the output from this table in the form of a vector of indicies which it can use
    output$datab = DT::renderDataTable(mrg(), server = F, rownames=F,
                                       options = list(searchHighlight = TRUE
                                                      ))
    #Some verbatim input for debugging purposes and to display to the user how many rows they have selected
    #Number of filtered rows is displayed in the table anyway
    output$tout = renderPrint({
        cat('\n\nFiltered Rows\n\n')
        #cat(input$datab_rows_all, sep = ', ')
        cat(length(input$datab_rows_all))
        cat('\n\nSelected rows:\n\n')
        #cat(input$datab_rows_selected, sep = ', ')
        cat(length(input$datab_rows_selected))
    })
    
    
    output$download <- downloadHandler(        
        filename <- function() {paste0("untitled.pdf")},
        content <- function(file) {
            pdf(file, width=input$dwidth, height=input$dheight)
            pMake()
            dev.off()
        }
    )
    
#     output$brushOut <- renderPrint({
#         brushedPoints(xy, yvar="gene", input$plot_brush)
#     })


    
}    
)
