library(shiny)
library(ggplot2)
library(biomaRt)


shinyServer(function(input,output,session){
  #funtions returning user inputs
  selectionX<-reactive({
    return(input$sampleX)
  })
  
  selectionY<-reactive({
    return(input$sampleY)
  })
  
  GOt<-reactive({
    return(input$GOterm)
  })
  
  productt<-reactive({
    return(input$productterm)
  })
  
  searchterm<-reactive({
    return(input$bygene)
  })
  
  WBIDs<-reactive({
    return(input$bywormgene)
  })
  
  dfxtype<-reactive({
    return(input$dfxtype)
  })
  
  dfytype<-reactive({
    return(input$dfytype)
  })
  
  inputnmin<-reactive({
    return(input$nMin)
  })
  
  

  
  output$distPlot <- renderPlot({
      #turning shiny user inputs into variables
      df1<-genewise_exp
      dafx<-dfxtype()
      dafy<-dfytype()
      selX<-selectionX()
      selY<-selectionY()
      dfx<-get(dafx)
      dfy<-get(dafy)
      Gts<-GOt()

      
      #base plot
      plot_out<-ggplot(df1,aes_string(x=dfx[selX],y=dfy[selY]))+geom_point()
      plot_out<-plot_out+labs(x=paste(dafx,selX),y=paste(dafy,selY))
      plot_out<-plot_out+theme_bw()
      
    
      
      #creating data frame for plotting given GOterm
      GO_data<-function(GOterm){
        genes<-get_genes_GO(GOterm)
        df<-data.frame(genes,selectionX=0,selectionY=0, GOterm,row.names=1)
        colnames(df)<-c(selX,selY,"ID")
        
        #with isoforms
        if(input$isoforms){
        for(gene in rownames(df)){
          for(i in rownames(dfx)){
            if(grepl(gene,i)){
                if(i %in% rownames(dfy)){
                  print(i)
                  df[i,1]<- dfx[i,selX]
                  df[i,2]<- dfy[i,selY]
                  df[i,3]<-GOterm
               }
            }
          }}}
        
        #without isoforms
        else{
        for(gene in rownames(df)){
          if ((gene %in% rownames(dfx))){
            if ((gene %in% rownames(dfy))){
            df[gene,1]<- dfx[gene,selX]
            df[gene,2]<- dfy[gene,selY]
          }
          }
        }}
        return(df)
      }

      #creating data frame for plotting given product term
      product_data<-function(keyterm){
        genes<-get_genes(keyterm,df_info,"product")
        df<-data.frame(genes,selectionX=0,selectionY=0, keyterm,row.names=1)
        colnames(df)<-c(selX,selY,"key")
        
        for(gene in rownames(df)){
          if ((gene %in% rownames(dfx))){
            if ((gene %in% rownames(dfy))){
              df[gene,1]<- dfx[gene,selX]
              df[gene,2]<- dfy[gene,selY]
            }
          }
        }
        return(df)
      }

 
    #create data frame for plotting WBIDs
      WBID_data<-function(WBIDs){
        genes<-getBM(attributes = c("refseq_mrna") , filters = c("ensembl_gene_id"), values=WBIDs, mart = ensembl)
        
        #1 means NR matches NR in our data
        for(i in 1:length(genes$refseq_mrna)){
          if(!(genes[i,"refseq_mrna"] %in% df_info$Name))
            genes[i,"NR"]<-0
          else
            genes[i,"NR"]<-1
        }
        
        NR_df<-subset(genes,genes$NR==1)
        df<-data.frame(NR_df$refseq_mrna,selectionX=0,selectionY=0,row.names=1)
        colnames(df)<-c(selX,selY)
        
        for(gene in rownames(df)){
          if ((gene %in% rownames(dfx))){
            if ((gene %in% rownames(dfy))){
              df[gene,1]<- dfx[gene,selX]
              df[gene,2]<- dfy[gene,selY]
            }
          }
        }
        return(df)
      }
      
      #create data frame with genes which aren't hightlighted
      not_high<-function(hlted,dfx,dfy){
        df<-data.frame(rownames(dfx),dfx[selX],dfy[selY],is_hi=1,row.names=1)
        for(i in rownames(hlted)){
          df[i,"is_hi"]<-0
        }
        df<-subset(df,df$is_hi==1)
        return(df)
      }
      
      
    
    #plot GO terms
     if(Gts!=""){
      GO_terms<-strsplit(Gts,",")[[1]]
      for(GO_ID in GO_terms){
        df2<-GO_data(GO_ID)
        plot_out<-plot_out+geom_point(data=df2, aes_string(x=df2[1],y=df2[2],colour="ID"))
      }
      }
      
    #plot product description key terms
    pterm<-productt()
      if(pterm!=""){
        p_terms<-strsplit(pterm,",")[[1]]
        for(pm in p_terms){
         df5<-product_data(pm)
         plot_out<-plot_out+geom_point(data=df5, aes_string(x=df5[1],y=df5[2],colour="key"))
        }
        }
    
    #highlight given gene yellow (or start of gene name)
    gene1<-searchterm()
    if(gene1!=""){
      genes1<-get_genes(gene1,df_info,"gene")
      for(gene in genes1){
          plot_out<-plot_out+geom_point(data=df1,aes_string(x=dfx[gene,selX],y=dfy[gene,selY]),colour="yellow", size=3.0)
        }
      
        }
  
    #highlight WBID red
    #plot product description key terms
    WBIDs<-WBIDs()
    if(WBIDs!=""){
      WBID_split<-strsplit(WBIDs," ")[[1]]
      df6<-WBID_data(WBID_split)
      View(df6)
      plot_out<-plot_out+geom_point(data=df6,aes_string(x=df6[1],y=df6[2]),colour="red")
      
      #not highlighted data frame
      df7<-not_high(df6,dfx,dfy)
      
      #p-value
      if(dafy=="genewise_tail_length" & dafx=="genewise_tail_length"){
        diff_else<-df7[selX]-df7[selY]
        diff_adele<-df6[selX]-df6[selY]
        p_val<-t.test(diff_else[1],diff_adele[1])
        p_val<-p_val[[3]] #get pvalue
        print(p_val)
        if(input$pval)
          plot_out<-plot_out+annotate(geom="text", x=20, y=90, label=paste0("p-value = ",round(p_val,digits=4)), color="black")
      }
    }
 

    
    #blue line
    if(dafx==dafy)
      plot_out<-plot_out+geom_abline(colour="blue",size=1.0)

    #make plot a square
    plot_out<-plot_out+coord_fixed()
    
    plot(plot_out)
    })
  
  #download plot
  output$downloadPlot<-downloadHandler(
    filename = function(){
      "file.pdf"
    },
    
    content = function(file){
      ggsave(file)
    }
  )
})





