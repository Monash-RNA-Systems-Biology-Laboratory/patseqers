library(shiny)
library(ggplot2)
library(biomaRt)


shinyServer(function(input,output,session){
  
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
  
  dataframx<-reactive({
    return(input$dataframx)
  })
  
  dataframy<-reactive({
    return(input$dataframy)
  })
  
  output$distPlot <- renderPlot({
    #ensembl = useMart("ensembl",dataset="scerevisiae_gene_ensembl")
    
      df1<-genewise_exp
      dafx<-dataframx()
      dafy<-dataframy()
      selX<-selectionX()
      selY<-selectionY()
      dfx<-get(dafx)
      dfy<-get(dafy)
      
      GO_data<-function(GOterm){
        genes<-get_genes_GO(GOterm)
        df<-data.frame(genes,selectionX=0,selectionY=0, GOterm,row.names=1)
        colnames(df)<-c(selX,selY,"ID")
        
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
          }}
        
#         for(gene in rownames(df)){
#           if ((gene %in% rownames(dfx))){
#             if ((gene %in% rownames(dfy))){
#             df[gene,1]<- dfx[gene,selX]
#             df[gene,2]<- dfy[gene,selY]
#           }
#           }
#         }
        return(df)
      }
      
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
      
      
    plot_out<-ggplot(df1,aes_string(x=dfx[selX],y=dfy[selY]))+geom_point()
    plot_out<-plot_out+labs(x=paste(dafx,selX),y=paste(dafy,selY))
    
    Gts<-GOt()
      
    #plot GO terms
     if(Gts!="None"){
      GO_terms<-strsplit(Gts,",")[[1]]
      for(GO_ID in GO_terms){
        df2<-GO_data(GO_ID)
        plot_out<-plot_out+geom_point(data=df2, aes_string(x=df2[1],y=df2[2],colour="ID"))
      }
      }
      
    #plot product description key terms
    pterm<-productt()
      if(pterm!="None"){
        p_terms<-strsplit(pterm,",")[[1]]
        for(pm in p_terms){
         df5<-product_data(pm)
         plot_out<-plot_out+geom_point(data=df5, aes_string(x=df5[1],y=df5[2],colour="key"))
        }
        }
    
      #one gene
    gene1<-searchterm()
    if(gene1!="Enter a gene"){
      genes1<-get_genes(gene1,df_info,"gene")
      for(gene in genes1){
        plot_out<-plot_out+geom_point(data=df1,aes_string(x=dfx[gene,selX],y=dfy[gene,selY]),colour="yellow", size=3.0)
      }
      }
    
      plot(plot_out)
    })
  
  output$downloadPlot<-downloadHandler(
    filename = function(){
      "file.pdf"
    },
    
    content = function(file){
      ggsave(file)
    }
      )
    
})





