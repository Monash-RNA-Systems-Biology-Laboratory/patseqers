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
  
  dfxtype<-reactive({
    return(input$dfxtype)
  })
  
  dfytype<-reactive({
    return(input$dfytype)
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

 
    #create data frame given file with genes to be hightlighted
      RBP_data<-function(){
        genes<-lapply(NR_df[,"refseq_mrna"],as.character)
        genes<-unlist(genes)
        product<-lapply(df_info[,"gene"],as.character)
        df<-data.frame(genes,selectionX=0,selectionY=0,row.names=1)
        colnames(df)<-c(selX,selY)
        
        for(i in ORF_df$wormbase_gene_seq_name){
          for(a in 1:length(product)){
            if(product[a] == i)
              genes<-c(genes,df_info[a,"Name"])}
        }
        
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
  
    #plot genes from file
    df6<-RBP_data()
    plot_out<-plot_out+geom_point(data=df6,aes_string(x=df6[1],y=df6[2]),colour="red")
   
    
    #p-value
    if(dafy=="genewise_tail_length" & dafx=="genewise_tail_length"){
      diff_all<-log2(dfx[selX]+0.5)-log2(dfy[selY]+0.5)
      diff_adele<-log2(df6[selX]+0.5)-log2(df6[selY]+0.5)
      p_val<-t.test(diff_all[1],diff_adele[1])
      p_val<-p_val[[3]] #get pvalue
      print(p_val)
      plot_out<-plot_out+annotate(geom="text", x=20, y=90, label=paste0("p-value = ",round(p_val,digits=4)), color="black")
      }
    
    #blue line
    if(dafx==dafy)
      plot_out<-plot_out+geom_abline(colour="blue",size=1.0)
    
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





