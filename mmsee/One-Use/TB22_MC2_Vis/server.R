#TB2 is shorter

library(reshape2)
library(ggplot2)

#Read in
readcountMC2 <- read.csv("readcountMC2.csv")
readcountTB2 <- read.csv("readcountTB2.csv")

#Set col_1 to rownames & Clip off row 1
rownames(readcountMC2) <- readcountMC2[,1]
rownames(readcountTB2) <- readcountTB2[,1]

readcountMC2[,1] <- NULL
readcountTB2[,1] <- NULL

shinyServer(function(input,output){
    
    #TB22_20
    TB2_WT <- data.frame(readcountTB2[,1:9])
    TB2_DS <- data.frame(readcountTB2[,10:18])
    #MC2_64_66
    MC2_WT <- data.frame(readcountMC2[,1:9])
    MC2_ST <- data.frame(readcountMC2[,10:18])
    MC2_SR <- data.frame(readcountMC2[,19:27])
    MC2_RR <- data.frame(readcountMC2[,28:36])
    
    TIME <- c(0,15,30,45,60,75,90,105,120)
    
    DF_LIST <- list(TB2_WT,TB2_DS,MC2_WT,MC2_ST,MC2_SR,MC2_RR)
    #View((DF_LIST))
    
    name.l <- rownames(readcountTB2)
    
    #tmp <- melt(t(TB2_WT))
    #tmp <- melt(t(readcountTB2))
    
    output$selgene <- renderUI({
        selectInput("genedisp", "Displaying gene", name.l)
    })
    
    output$normdisp <- renderUI({
        selectInput("normdisp", "Normalising against", c("SRP68", "PRP46"))
    })
    
    #Plot is normalised against genes
    mkplot1 <- reactive({
        #TB22
        tmp <- melt(t(TB2_WT))
        tmp$strain <- "TB22_WT"
        tmp$TIME <- TIME
        normg <- tmp[tmp$Var2 == input$normdisp,]
        tmp <- tmp[tmp$Var2 == input$genedisp,]
        tmp$strnorm <- (tmp$value/normg$value)*mean(normg$value)
        tmp$strnorm <- tmp$strnorm/mean(tmp$strnorm)
        
        tmp2 <- melt(t(TB2_DS))
        tmp2$strain <- "TB22_DS"
        tmp2$TIME <- TIME
        normg <- tmp2[tmp2$Var2 == input$normdisp,]
        tmp2 <- tmp2[tmp2$Var2 == input$genedisp,]
        tmp2$strnorm <- (tmp2$value/normg$value)*mean(normg$value)
        tmp2$strnorm <- tmp2$strnorm/mean(tmp2$strnorm)
        
        #MC2_64_66
        tmp3 <- melt(t(MC2_WT))
        tmp3$strain <- "MC2_WT"
        tmp3$TIME <- TIME
        normg <- tmp3[tmp3$Var2 == input$normdisp,]
        tmp3 <- tmp3[tmp3$Var2 == input$genedisp,]
        tmp3$strnorm <- (tmp3$value/normg$value)*mean(normg$value)
        tmp3$strnorm <- tmp3$strnorm/mean(tmp3$strnorm)
        
        tmp4 <- melt(t(MC2_ST))
        tmp4$strain <- "MC2_SET1"
        tmp4$TIME <- TIME
        normg <- tmp4[tmp4$Var2 == input$normdisp,]
        tmp4 <- tmp4[tmp4$Var2 == input$genedisp,]
        tmp4$strnorm <- (tmp4$value/normg$value)*mean(normg$value)
        tmp4$strnorm <- tmp4$strnorm/mean(tmp4$strnorm)
        
        tmp5 <- melt(t(MC2_SR))
        tmp5$strain <- "MC2_SET1_RRP6"
        tmp5$TIME <- TIME
        normg <- tmp5[tmp5$Var2 == input$normdisp,]
        tmp5 <- tmp5[tmp5$Var2 == input$genedisp,]
        tmp5$strnorm <- (tmp5$value/normg$value)*mean(normg$value)
        tmp5$strnorm <- tmp5$strnorm/mean(tmp5$strnorm)
        
        tmp6 <- melt(t(MC2_RR))
        tmp6$strain <- "MC2_RRP6"
        tmp6$TIME <- TIME
        normg <- tmp6[tmp6$Var2 == input$normdisp,]
        tmp6 <- tmp6[tmp6$Var2 == input$genedisp,]
        tmp6$strnorm <- (tmp6$value/normg$value)*mean(normg$value)
        tmp6$strnorm <- tmp6$strnorm/mean(tmp6$strnorm)
        
        plot <- ggplot(data=tmp, aes(x=TIME, y=strnorm)) + geom_line(aes(color=strain))
        plot <- plot + geom_line(data=tmp2, aes(x=TIME, y=strnorm, color=strain))
        plot <- plot + geom_line(data=tmp3, aes(x=TIME, y=strnorm, color=strain))
        plot <- plot + geom_line(data=tmp4, aes(x=TIME, y=strnorm, color=strain))
        plot <- plot + geom_line(data=tmp5, aes(x=TIME, y=strnorm, color=strain))
        plot <- plot + geom_line(data=tmp6, aes(x=TIME, y=strnorm, color=strain))
        #Plot Style
        plot <- plot + scale_x_continuous(breaks=TIME) + scale_y_continuous(breaks=c(0:100))
        plot <- plot + theme_bw() + theme(
            text=element_text(size=15),
            legend.justification=c(1,1), 
            legend.position=c(1,1))
        plot <- plot + labs(x="Time (minutes)", y="Scaled counts")
        return(plot)
    })
    output$plot <- renderPlot({
            mkplot1()
    })
    output$PlotD1 <- downloadHandler(
        filename = function() { paste0('Untitled.pdf', sep='') },
        content = function(file) {
            ggsave(file,mkplot1())
        })
})