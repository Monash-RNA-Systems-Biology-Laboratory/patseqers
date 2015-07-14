

#server.R
#Authors: Michael See, Andrew Pattinson

library(ggplot2)
library(edgeR) 
library(reshape2)
library(plyr)


shinyServer(function(input,output) {
    
    rdcsv <- reactive({
        orig.df <- data.frame(read.csv(paste0(input$datset,"-count.csv"), header = TRUE))
        name.df <- data.frame(read.csv(paste0(input$datset,"-info.csv"), header = TRUE))
        
        nameN.df <- data.frame(name.df[,3], name.df[,1])
        colnames(nameN.df) <- c("Gene", "Name")
        
        origN.df <- merge(orig.df, nameN.df)
        nameN.df <- merge(orig.df, nameN.df)
        #origN.df <- subset(origN.df, !duplicated(origN.df$Gene))
        dup <- duplicated(origN.df$Gene)
        nameN.df <- subset(nameN.df, !duplicated(origN.df$Gene))
        orig.df <- subset(orig.df, !duplicated(origN.df$Gene))
        
        origN.df <- orig.df[,-1]
        rownames(origN.df) <- nameN.df$Gene #Gene names are now appropriately added as rownames
        
        #Load in config file
        
        return(origN.df)
    })
    
    rdcfg <- reactive({
        return(cfg.df <- read.csv(paste0(input$datset,"-confg.csv")))
    })
    
    
    numGenGraph <- reactive({
        return(length(input$geneToGraph))
    })
    minVar <- reactive({
        return(input$nMin)
    })
    gtn <- reactive({
        return(input$geneToNorm)
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
                data.frame(rdcsv()[apply(rdcsv(),1,function(row) {
                    any(row >= input$nMin)
                }),])
        })
    
    nPlot <- reactive({
        input$nPlotV
    })
    
    output$selDataSet <- renderUI({
        sets.df <- read.csv("setname.csv", header=T)
        sets.df <- data.frame(lapply(sets.df, as.character), stringsAsFactors=F)
        stvec <- as.character(sets.df$setname)
        selectInput("datset", "Choose dataset",sets.df$setname)
    })
    
    output$SelGeneNorm <- renderUI({
        name.l <- rownames(flt())
        selectInput("geneToNorm", "Choose gene to normalise against", name.l)
    })
    
    output$selGene <- renderUI({
        name.l <- rownames(flt())
        selectizeInput("geneToGraph", "Choose a gene to graph", choices = name.l, multiple = F, selected=name.l[2])
    })
    
    observeEvent(input$do,{
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
            
            genetograph.ls <- rep(list(data.frame(0)),len)
            for(i in 1:len){
                genetograph.ls[[i]] <- subset(mlt, Name==paste0(namelist[i]))
            }
            
            graphlist.ls <- rep(list(list(0)),len)
            for(i in 1:len){
                graphlist.ls[[i]] <- rdcfg()
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
            
            return(graphlist.ls)
        }
        
        
        
        tmmFunc <- function(indf){
            flt.df <- indf
            flt.df <- DGEList(flt.df)
            flt.df <- calcNormFactors(flt.df)
            ret <- data.frame(cpm(flt.df, log=T, prior.count=input$pRC))
            return(ret)
        }
        tmmFunc2 <- function(indf){
            flt.df <- indf
            flt.df <- DGEList(flt.df)
            flt.df <- calcNormFactors(flt.df)
            ret <- data.frame(cpm(flt.df, log=F, prior.count=input$pRC))
            return(ret)
        }
        
        plotHelper2 <- function(m2){
            plot2 <- ggplot()
            for(i in 1:numGenGraph()){
                for(j in 1:length(m2[[1]])){
                    m2[[i]][[j]]$linename <- (merge(m2[[i]][[j]], rdcfg(), "line")$linename[1])
                    #m2[[i]][[j]]$lineName <- paste0(m2[[i]][[j]]$Name, "_Line",m2[[i]][[j]]$line)
                    m2[[i]][[j]]$unnormVarl2 <- log2(m2[[i]][[j]]$value)
                    plot2 <- plot2 + geom_line(data=m2[[i]][[j]], aes(x=xaxis,y=unnormVarl2, colour=linename))
                }
            }
            return(plot2)
        }
        
        fldcng <- function(df){
            a <- which(df$fldcng!=0, arr.ind=T)
            df$foldvalue <- df$value/df$value[a[1]]
            for(i in 1:length(df$foldvalue)){
                if(df$foldvalue[i] < 1){
                    df$foldvalue[i] <- 1/df$foldvalue[i]*-1 
                }
            }
            return(df)
        }
        
        fldcng2 <- function(df){
            a <- which(df$fldcng!=0, arr.ind=T)
            df$foldvaluetmm <- df$normVal/df$normVal[a[1]]
            for(i in 1:length(df$foldvaluetmm)){
                if(df$foldvaluetmm[i] < 1){
                    df$foldvaluetmm[i] <- 1/df$foldvaluetmm[i]*-1 
                }
            }
            return(df)
        }
        
        plotHelper2siz1 <- function(m2){
            #m2 data to graph
            
            tp.df <- data.frame(m2[[1]][[1]])
            for(i in 1:length(m2)){
                for(j in 1:length(m2[[1]])){
                    tp.df <- rbind(tp.df, m2[[i]][[j]])
                }
            }
            tp.df <- tp.df[-1,]
            rd <- rdcfg()
            rd <- rd[,-3:-5]
            tp.df$dub <- floor(tp.df$X)
            g1 <- merge(rd, tp.df, by.x="X", by.y="dub")
            g1 <- fldcng(g1)

            plot2 <- ggplot(g1, aes(x=linename,y=foldvalue, fill=factor(linename))) + geom_bar(stat="identity", width = 0.75)
            plot2 <- plot2 + ggtitle(paste0("Fold change in relation to ", g1$linename[which(g1$fldcng!=0, arr.ind=T)]))
            
            return(plot2)
        }
        
        
        plotInput2 <- reactive({
            tograph <- meltsplit(flt(), gtg())
            
            if(dim(tograph[[1]][[1]]) == 1){
                p2 <- plotHelper2siz1(tograph)
                p2 <- p2 + ylab("Fold change") + xlab("Strain")
            } else {
                p2 <- plotHelper2(meltsplit(flt(), gtg()))
                p2 <- p2 + ggtitle(paste0("Log2 of un-normalised counts")) + ylab("Log2(RPM)") + xlab("Time")
            }
            return(p2)
        })
        
        
        output$plot2 <- renderPlot({
            return(plotInput2())
        })
        
        plotHelper3 <- function(m3,m2){
            plot3 <- ggplot()
            
            tn  <- meltsplit(flt(), gtn())
            
            if(input$drTMM){
                for(i in 1:numGenGraph()){
                    for(j in 1:length(m3[[1]])){
                        m3[[i]][[j]]$linename <- (merge(m3[[i]][[j]], rdcfg(), "line")$linename[1])
                        #m3[[i]][[j]]$lineName <- paste0(m3[[i]][[j]]$Name, "_Line",m3[[i]][[j]]$line)
                        m3[[i]][[j]]$normVarl2 <- m2[[1]][[j]]$value
                        plot3 <- plot3 + geom_line(data=m3[[i]][[j]], aes(x = xaxis, y = normVarl2, colour=linename))
                    }
                }
            } else {
                for(i in 1:numGenGraph()){
                    for(j in 1:length(m3[[1]])){
                        m3[[i]][[j]]$linename <- (merge(m3[[i]][[j]], rdcfg(), "line")$linename[1])
                        #m3[[i]][[j]]$lineName <- paste0(m3[[i]][[j]]$Name, "_Line",m3[[i]][[j]]$line)
                        m3[[i]][[j]]$naive <- m3[[i]][[j]]$value/tn[[1]][[j]]$value
                        plot3 <- plot3 + geom_line(data=m3[[i]][[j]], aes(x = xaxis, y = naive, colour=linename)) +
                            ylab("log2 RPM") + xlab("Time")
                    }
                }
            }
            return(plot3)
        }
        
        plotHelper3siz1 <- function(m3,m2){
            #m3 data to graph, m2 normalisation data
            
            tp.df <- data.frame(m3[[1]][[1]])
            tn.df <- data.frame(m2[[1]][[1]])
            for(i in 1:length(m3)){
                for(j in 1:length(m3[[1]])){
                    tp.df <- rbind(tp.df, m3[[i]][[j]])
                }
            }
            
            for(i in 1:length(m2)){
                for(j in 1:length(m2[[1]])){
                    tn.df <- rbind(tn.df, m2[[i]][[j]])
                }
            }
            
            nrn <- meltsplit(flt(), gtn())
            
            tmp.df <- data.frame(m2[[1]][[1]])
            for(i in 1:length(nrn)){
                for(j in 1:length(nrn[[1]])){
                    tmp.df <- rbind(tmp.df, nrn[[i]][[j]])
                }
            }
            
            tp.df <- tp.df[-1,]
            tn.df <- tn.df[-1,]
            tmp.df <- tmp.df[-1,]
            tmp.df <- tmp.df[-1*dim(tmp.df)[1],]
            tp.df$normVal <- tn.df$value
            tp.df$value <- tp.df$value/tmp.df$value
            rd <- rdcfg()
            rd <- rd[,-3:-5]
            tp.df$dub <- floor(tp.df$X)
            g1 <- merge(rd, tp.df, by.x="X", by.y="dub")
            g1 <- fldcng(g1)
            g1 <- fldcng2(g1)
           
            if(input$drTMM){
                plot3 <- ggplot(g1, aes(x=linename,y=foldvaluetmm, fill=factor(linename))) + geom_bar(stat="identity", width = 0.75)
                plot3 <- plot3 + ggtitle(paste0("Fold change in relation to ", g1$linename[which(g1$fldcng!=0, arr.ind=T)], " After TMM normalisation"))
            } else {
                plot3 <- ggplot(g1, aes(x=linename,y=foldvalue, fill=factor(linename))) + geom_bar(stat="identity", width = 0.75)
                plot3 <- plot3 + ggtitle(paste0("Fold change in relation to ", g1$linename[which(g1$fldcng!=0, arr.ind=T)], " Normalised against ", gtn()))
            }
            
            return(plot3)
        }
        
        plotInput3 <- reactive({
            tograph <- meltsplit(flt(), gtg())
            
            if(dim(tograph[[1]][[1]]) == 1){
                p3 <- plotHelper3siz1(tograph ,meltsplit(tmmFunc2(flt()), gtg()))
                if(input$drTMM){
                    p3 <- p3 + ylab("TMM (Fold Change)") + xlab("Strain")
                } else {
                    p3 <- p3 + ylab(paste0("Fold Change Normalised against single gene (", gtn(), ")")) + xlab("Strain")
                }
            } else {
                p3 <- plotHelper3(tograph ,meltsplit(tmmFunc(flt()), gtg()))
                if(input$drTMM){
                    p3 <- p3 + ggtitle(paste0(gtg(),"Against library")) + ylab("TMM (RPM)") + xlab("Time")
                } else {
                    p3 <- p3 + ggtitle(paste0(gtg(),"/ Against single gene (", gtn(), ")")) + ylab("Log2 RPM") + xlab("Time")
                }
            }
            
            return(p3)
        })
        output$plot3 <- renderPlot({
            return(plotInput3())
        })
        
        
        output$PlotD1 <- downloadHandler(
            filename = function() { paste0(input$dataset, 'Untitled.pdf', sep='') },
            content = function(file) {
                ggsave(file,plotInput2())
            })
        
        
        output$PlotD2 <- downloadHandler(
            filename = function() { paste0(input$dataset, 'Untitled.pdf', sep='') },
            content = function(file) {
                ggsave(file,plotInput3())
        
            })
        
    })
})


