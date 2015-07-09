#server.R
#Authors: Michael See, Andrew Pattinson

#Imports
library(edgeR)
library(ggplot2)

library(reshape2)
library(plyr)

###Load in the .csv files to be plotted. orig.df is the counts, name.df is just for grabbing the gene names###
orig.df <- data.frame(read.csv("genewise-count.csv", header = TRUE))
name.df <- data.frame(read.csv("genewise-info.csv", header = TRUE))

nameN.df <- data.frame(name.df[,3], name.df[,1])
colnames(nameN.df) <- c("Gene", "Name")

origN.df <- merge(orig.df, nameN.df)
nameN.df <- merge(orig.df, nameN.df)

origN.df <- orig.df[,-1]


rownames(origN.df) <- nameN.df$Gene #Gene names are now appropriately added as rownames
#View(origN.df)
#Load in config file
cfg.df <- read.csv("confg.csv")
#View(cfg.df)

###END LOADING CSV ###
#View(origN.df)
shinyServer(function(input,output) {
    
    numGenGraph <- reactive({
        return(length(input$geneToGraph))
    })
    minVar <- reactive({
        input$nMin
    })
    gtn <- reactive({
        input$geneToNorm
    })
    
    output$n <- reactive({
        return(input$nPlotV)
    })
    gtg <- reactive({
        return(input$geneToGraph)
        
    })

    flt <-
        reactive({
            flt.df <-
                data.frame(origN.df[apply(origN.df,1,function(row) {
                    any(row >= input$nMin)
                }),])
        })
    
    nPlot <- reactive({
        input$nPlotV
    })
    
    output$SelGeneNorm <- renderUI({
        name.l <- rownames(flt())
        selectInput("geneToNorm", "Choose gene to normalise against", name.l)
    })
    
    output$selGene <- renderUI({
        name.l <- rownames(flt())
        selectizeInput("geneToGraph", "Choose a gene to graph", choices = name.l, multiple = T, selected=name.l[2])
    })
    
    observeEvent(input$do,{
        View(flt())
        if(is.null(gtg())){
            print(gtg())
            return(0)
        }
        
        #numDup <- input$numRep + 1 #number of duplicates done
        #numSet <- input$setNum #Number of sets of data
        #timVec <- unique(na.omit(as.numeric(unlist(strsplit(unlist(input$timScale), "[^0-9]+"))))) Depricated
        
        meltsplit <- function(tomelt, namelist){
            tomelt$Name <- rownames(tomelt)
            mlt <- melt(tomelt, "Name")
            len <- length(namelist)
            #View(mlt)
            
            genetograph.ls <- rep(list(data.frame(0)),len)
            for(i in 1:len){
                genetograph.ls[[i]] <- subset(mlt, Name==paste0(namelist[i]))
            }
            
            graphlist.ls <- rep(list(list(0)),len)
            for(i in 1:len){
                graphlist.ls[[i]] <- cfg.df
            }
            #rep(list(list(0)),len)
            
            for(i in 1:len){
                graphlist.ls[[i]] <- merge(genetograph.ls[[i]], graphlist.ls[[i]], by.x="variable", by.y="Sample")
            }
            for(i in 1:len){
                graphlist.ls[[i]] <- split(graphlist.ls[[i]], graphlist.ls[[i]]$line)
            }
            for(i in 1:length(graphlist.ls)){
                #print(paste0(len,"a"))
                
                for(j in 1:length(graphlist.ls[[i]])){
                    graphlist.ls[[i]][[j]] <- ddply(graphlist.ls[[i]][[j]], "xaxis", numcolwise(mean))
                    graphlist.ls[[i]][[j]]$Name <- genetograph.ls[[i]]$Name[1]
                }
                
            }
            #View(graphlist.ls)
            
            return(graphlist.ls)
        }
        

        
        tmmFunc <- function(indf){
            flt.df <- indf
            flt.df <- DGEList(flt.df)
            flt.df <- calcNormFactors(flt.df)
            ret <- data.frame(cpm(flt.df, log=T, prior.count=input$pRC))
            return(ret)
        }

        
        plotHelper2 <- function(m2){
            plot2 <- ggplot()
            for(i in 1:numGenGraph()){
                for(j in 1:length(m2[[1]])){
                    m2[[i]][[j]]$lineName <- paste0(m2[[i]][[j]]$Name, "_Line",m2[[i]][[j]]$line)
                    m2[[i]][[j]]$unnormVarl2 <- log2(m2[[i]][[j]]$value)
                    plot2 <- plot2 + geom_line(data=m2[[i]][[j]], aes(x=xaxis,y=unnormVarl2, colour=lineName))
                }
            }
            return(plot2)
        }
        plotInput2 <- reactive({
            p2 <- plotHelper2(meltsplit(flt(), gtg()))
            p2 <- p2 +ggtitle(paste0("Log2 of unormalised counts")) + ylab("Log2(RPM)") + xlab("Time(minutes)")
            return(p2)
        })
        output$plot2 <- renderPlot({
            return(plotInput2())
        })
        
        plotHelper3 <- function(m3,m2){
            plot3 <- ggplot()
            for(i in 1:numGenGraph()){
                for(j in 1:length(m3[[1]])){
                    m3[[i]][[j]]$lineName <- paste0(m3[[i]][[j]]$Name, "_Line",m3[[i]][[j]]$line)
                    m3[[i]][[j]]$normVarl2 <- log2((m3[[i]][[j]]$value / m2[[1]][[j]]$value))
                    plot3 <- plot3 + geom_line(data=m3[[i]][[j]], aes(x=xaxis,y=normVarl2, colour=lineName))
                }
            }
            return(plot3)
        }
        plotInput3 <- reactive({
            p3 <- plotHelper3(meltsplit(flt(), gtg()),meltsplit(tmmFunc(flt()), gtn()))
            p3 <- p3 + ggtitle(paste0(gtg(),"/ TMM normalisation")) + ylab("TMM (RPM)") + xlab("Time(minutes)")
            return(p3)
        })
        output$plot3 <- renderPlot({
            return(plotInput3())
        })
        
    })
})