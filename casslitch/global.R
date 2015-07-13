raw_gene_exp<-read.csv("genewise-count.csv")
df_rawgene_exp<-data.frame(raw_gene_exp,  stringsAsFactors = T)
df3<-df_rawgene_exp[,-1]
rownames(df3)<-df_rawgene_exp[,1]

gene_nwr<-colnames(df3)
gene_n<-list()
for(i in gene_nwr){
  gen<-substr(i,0,nchar(i)-5)
  if((gen %in% gene_n)==FALSE)
    gene_n[length(gene_n)+1]<-gen
}

for(a in gene_n){
  var1<-paste0(a,".rep1")
  print(var1)
  var2<-paste0(a,".rep2")
  m<-paste0(a,".mean")
  df3[,m]<-((df3[,var1]+df3[,var2])/2)
}


df.log<-log2(df3[,]+1)
df1<-df.log
check_boxes<-names(df1)