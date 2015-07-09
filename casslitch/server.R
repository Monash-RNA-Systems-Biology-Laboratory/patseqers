library(shiny)
library(ggplot2)
library(biomaRt)

# raw_gene_exp<-read.csv("genewise-count.csv")
# df_rawgene_exp<-data.frame(raw_gene_exp,  stringsAsFactors = T)
# df3<-df_rawgene_exp[,-1]
# rownames(df3)<-df_rawgene_exp[,1]

# gene_nwr<-colnames(df1)
# gene_n<-list()
# for(i in gene_nwr){
#   gen<-substr(i,0,nchar(i)-5)
#   if((gen %in% gene_n)==FALSE)
#     gene_n[length(gene_n)+1]<-gen
# }
# 
# for(a in gene_n){
#   var1<-paste0(a,".rep1")
#   print(var1)
#   var2<-paste0(a,".rep2")
#   m<-paste(a,"mean")
#   df1[,m]<-((df1[,var1]+df1[,var2])/2)
# }


# df.log<-log2(df3[,]+1)
# df1<-df.log

shinyServer(function(input, output,session){
    selectionX<-reactive({
      return(input$sampleX)
    })
    
    selectionY<-reactive({
      return(input$sampleY)
    })
    
   GOt<-reactive({
      return(input$GOterm)
    })
  
  searchterm<-reactive({
    return(input$bygene)
  })
    
  output$distPlot <- renderPlot({
      ensembl = useMart("ensembl",dataset="scerevisiae_gene_ensembl")
      selX<-selectionX()
      selY<-selectionY()
      
#       if(grepl("mean",selX))
#         selX<-substr(selX,1,nchar(selX)-5)
#   
#       if(grepl("mean",selY))
#         selY<-substr(selY,1,nchar(selY)-5)
#       
  
      GO_data<-function(GOterm){
        GO_t<-GOterm
        genes<-getBM(attributes = c("ensembl_gene_id") , filters = c("go_id"), values = "GO:0005634" , mart = ensembl)
        df<-data.frame(genes,selectionX=0,selectionY=0, GO_t,row.names=1)
        colnames(df)<-c(selX,selY,"ID")
        
        for(gene in rownames(df)){
          if ((gene %in% rownames(df1))){
            df[gene,selX]<- df1[gene,selX]
            df[gene,selY]<- df1[gene,selY]
          }
        }
        return(df)
      }
      
      
    plot_out<-ggplot(df1,aes_string(x=selX,y=selY))+geom_point()
    
    Gts<-GOt()
      
     if(Gts!="None"){
      GO_terms<-strsplit(Gts,",")
      for(GO_ID in GO_terms){
        df2<-GO_data(GO_ID)
        plot_out<-plot_out+geom_point(data=df2, aes_string(x=selX,y=selY,colour="ID"))
      }
      }
      
      #one gene
      gene1<-searchterm()
      data1<-df1[gene1,]
      plot_out<-plot_out+geom_point(data=data1,aes_string(x=selX,y=selY),colour="yellow", size=3.0)
      
      
      plot(plot_out)
    })
    
})