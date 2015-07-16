library(shiny)
library(ggplot2)
library(biomaRt)


shinyServer(function(input,output,session){
  #make inputs reactive
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
  
  #Plots
  output$distPlot <- renderPlot({
    
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
      
      else{
        #without isoforms       
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
    
    #xmax
    xmax<-max(dfx[selX],na.rm = T)
    #ymax
    ymax<-max(dfy[selY],na.rm = T)
    
    plot_out<-ggplot(df1,aes_string(x=dfx[selX],y=dfy[selY]))+geom_point()
    plot_out<-plot_out+labs(x=paste(dafx,selX),y=paste(dafy,selY))
    plot_out<-plot_out+geom_segment(aes_string(x=0,y=0,xend=xmax,yend=ymax),colour="blue",size=1.7)
    
    Gts<-GOt()
    
    #plot GO term(s)
    if(Gts!=""){
      GO_terms<-strsplit(Gts,",")[[1]]
      for(GO_ID in GO_terms){
        df2<-GO_data(GO_ID)
        plot_out<-plot_out+geom_point(data=df2, aes_string(x=df2[1],y=df2[2],colour="ID"))
      }
    }
    
    #plot by product description key term(s)
    pterm<-productt()
    if(pterm!=""){
      p_terms<-strsplit(pterm,",")[[1]]
      for(pm in p_terms){
        df5<-product_data(pm)
        plot_out<-plot_out+geom_point(data=df5, aes_string(x=df5[1],y=df5[2],colour="key"))
      }
    }
    
    #plot by gene/gene family
    gene1<-searchterm()
    if(gene1!=""){
      genes1<-get_genes(gene1,df_info,"gene")
      dfz<-data.frame(genes1,xax=0,yax=0,row.names=genes1)
      for(gene in genes1){
        dfz[gene,"xax"]<-dfx[gene,selX]
        dfz[gene,"yax"]<-dfy[gene,selY]
      }
      plot_out<-plot_out+geom_point(data=dfz,aes_string(x="xax",y="yax"),colour="yellow", size=3.0)
      
    }
    
    
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





