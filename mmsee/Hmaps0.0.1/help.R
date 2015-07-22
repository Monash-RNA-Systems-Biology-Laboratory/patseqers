#help.R
#Authors - Paul Harrison
#        - Michael See

t1 <- function(filename.prefix, elist, 
         min.sd=0.0, min.span=0.0, min.svd=0.0, svd.rank=NULL,
         annotation1=c('gene'),
         annotation2=c('product'), 
         res=150, row.labels=NA,
         reorder.columns = FALSE) {    
    # Don't die if missing annotation
    annotation1 <- annotation1[ annotation1 %in% colnames(elist$gene) ]
    annotation2 <- annotation2[ annotation2 %in% colnames(elist$gene) ]
    keep <- rep(TRUE, nrow(elist$E))
    
    span <- row.apply(elist$E, max) - row.apply(elist$E, min)
    keep <- (keep & span > 0.0) #Always filter genes that are exactly identical in all
    if (min.span > 0.0) {
        keep <- (keep & span >= min.span)
    }
    
    if (min.sd > 0.0) {
        sd <- sqrt(row.apply(elist$E, var))
        keep <- (keep & sd >= min.sd)
    }
    
    if (min.svd > 0.0) {
        #if (is.null(svd.rank))
        #    svd.rank <- ncol(elist$E)-1
        #s <- svd(t(scale(t(elist$E), center=TRUE,scale=FALSE)), nu=svd.rank,nv=svd.rank)
        #cat('SVD d diagonal:\n')
        #print(s$d[basic.seq(svd.rank)])
        #mag <- sqrt( rowSums(s$u*s$u) * nrow(s$u) / ncol(s$u) )
        #keep <- (keep & mag >= min.svd)
        
        keep <- keep & svd.gene.picker(elist$E, svd.rank, min.svd)
    }
    
    elist <- elist[keep,]
    
    #cat(sprintf("%d genes shown\n", nrow(elist)))
    
    averages <- rowMeans(elist$E)
    
    if (is.na(row.labels))
        row.labels <- (nrow(elist) <= 600)
    
    data <- t(scale(t(elist$E), center=TRUE,scale=FALSE))
    
    labels <- list()
    
    for(colname in annotation1)
        if (!all(is.na(elist$gene[,colname])))
            labels[[ length(labels)+1 ]] <- elist$gene[,colname]
    
    labels[[ length(labels)+1 ]] <- rownames(data)
    
    for(colname in annotation2)
        if (!all(is.na(elist$gene[,colname])))
            labels[[ length(labels)+1 ]] <- elist$gene[,colname]
    
    height <- if(row.labels) (10*nrow(data)+1500)*res/150 else 2500*res/150
    #png(sprintf('%s.png',filename.prefix), width=2000*res/150, height=height, res=res)
    
    #heatmap <- nesoni.heatmap(data, col=signed.col, symkey=TRUE,symbreaks=TRUE, labRow=(if(row.labels) NULL else NA), margins=margins, main=main, ...)
    heatmap <- nesoni.heatmap(data, labels=if(row.labels) labels else list(), reorder.columns=reorder.columns, levels=averages)
        #dev.off()

    shuffled.elist <- elist[rev(heatmap$dend.row$order),]
    
    table.filename <- sprintf('%s.csv', filename.prefix)
    
    sink(table.filename)
    cat('# Heatmap data\n')
    cat('#\n')
    cat('# Values given are log2 reads per million\n')
    cat('#\n')
    cat(sprintf('# %d genes shown\n', nrow(data)))
    cat('#\n')
    
    frame <- data.frame(name=rownames(shuffled.elist$E), row.names=rownames(shuffled.elist$E), check.names=FALSE) 
    for(colname in c(annotation1,annotation2)) {
        frame[,colname] <- shuffled.elist$gene[,colname]
    }     
    #frame[,'cluster hierarchy'] <- rev(dendrogram.paths(heatmap$rowDendrogram))
    frame[,'cluster hierarchy'] <- rev(heatmap$dend.row$paths)
    frame <- data.frame(frame, shuffled.elist$E, check.names=FALSE)
    
    write.csv(frame, row.names=FALSE)
    sink()
    
    invisible(list(heatmap=heatmap, frame=frame))
}