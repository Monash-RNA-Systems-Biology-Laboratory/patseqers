#Author: Cassandra Litchfield
#clitchfield@student.unimelb.edu.au
#July 2015

library(biomaRt)
ensembl = useMart("ensembl",dataset="scerevisiae_gene_ensembl")

#csvs
df_info<-read.csv("genewise-info.csv")
raw_gene_exp<-read.csv("genewise-count.csv")
raw_gene_tail<-read.csv("genewise-tail.csv")

#delete duplicates
df_info<-subset(df_info,!duplicated(df_info$gene))
raw_gene_exp<-subset(raw_gene_exp,!duplicated(raw_gene_exp$Name))
raw_gene_tail<-subset(raw_gene_tail,!duplicated(raw_gene_tail$Name))

#df info
rownames(df_info)<-df_info[,3]

#gene expression
df_rawgene_exp<-data.frame(raw_gene_exp,  stringsAsFactors = T)
df3<-df_rawgene_exp[,-1]
rownames(df3)<-df_rawgene_exp[,1]

#get sample names
samples<-colnames(df3)
for(i in 1:length(samples))
  samples[i]<-strsplit(samples[i],split=".rep")[[1]]
samples<-subset(samples,!duplicated(samples))

#count number of replicates for given sample
sample_replicates<-data.frame(samples,replicates=0,row.names = samples)
for(i in samples){
  reps<-subset(colnames(df3),grepl(paste0(i,".rep"),colnames(df3)))
  sample_replicates[i,"replicates"]<-length(reps)
}

#calculate means for each sample
original_names<-names(df3)
for(i in samples){
  reps<-subset(colnames(df3),grepl(paste0(i,".rep"),colnames(df3)))
  df8<-subset(df3, select=reps)
  df3<-data.frame(df3,rowMeans(df8))
}

col_nom<-c(original_names,lapply(samples,paste0,"_mean"))
colnames(df3)<-col_nom


#take log2 of all coloumns
genewise_exp<-log2(df3[,]+1)


#tail lengths
df_rawgene_tail<-data.frame(raw_gene_tail,  stringsAsFactors = T)
df4<-df_rawgene_tail[,-1]
rownames(df4)<-df_rawgene_tail[,1]

#calculate means for each sample (tail lengths)
for(i in samples){
  reps<-subset(colnames(df3),grepl(paste0(i,".rep"),colnames(df3)))
  df8<-subset(df4, select=reps)
  df4<-data.frame(df4,rowMeans(df8))
}

colnames(df4)<-col_nom

genewise_tail_length<-df4


#drop down boxes
check_boxes<-names(df3)

#returns genes with given product description key term
get_genes<-function(key_term,dfinfo,colu) {
  genes_bool<-rep(F,nrow(dfinfo))
  for(i in 1:nrow(dfinfo)){
    if(grepl(key_term, dfinfo[i,colu])[1])
      genes_bool[i]<-T
  }
  genes<-rep("",sum(genes_bool, na.rm=TRUE))
  count=1
  for(i in 1:length(genes_bool)){
    if(genes_bool[i]){
      genes[count]<-as.character(dfinfo[i,"Name"])
      count<-count+1
    }
  }
  print(length(genes))
  return(genes)
}

#returns genes for given GO_term
get_genes_GO<-function(GO_term){
  genes<-getBM(attributes = c("ensembl_gene_id") , filters = c("go_id"), values = GO_term , mart = ensembl)
  return(genes)
}