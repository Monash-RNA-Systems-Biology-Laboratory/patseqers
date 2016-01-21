setwd("/data/home/apattison/Bioinformatics/mouse_breast_cancer/feature_counts/sam_and_feature_files")

modify_gff_inplace <- function (gff_file) {
  start_gff_file <- read.delim(gff_file, header=FALSE,
                               comment.char="",stringsAsFactors=F)
  
  colnames(start_gff_file)<- c('Chromosome', 'Generated_By', 'Feature_Type', 
                               'Peak_Start','Peak_End','-',
                               'Orientation', '--','Information')
  
  plus_frame <- start_gff_file[start_gff_file[,'Orientation'] == '+',]
  plus_frame [,c('Peak_Start', 'Peak_End')] <- plus_frame [,c('Peak_Start', 'Peak_End')]+12
  
  plus_reads <- plus_frame[
    with(plus_frame,order(
      Chromosome,Orientation,Peak_Start)
    ),
    ]
  
  minus_frame <- start_gff_file[start_gff_file[,'Orientation'] == '-',]
  minus_frame [,c('Peak_Start', 'Peak_End')] <- minus_frame [,c('Peak_Start', 'Peak_End')]-12
  
  minus_reads<- minus_frame[
    with(minus_frame,order(
      Chromosome,Peak_Start)
    ),
    ]
  for (row in 1:nrow(plus_reads)){
    if (row == 1){
      next
    }
    if (plus_reads[row, 'Chromosome'] != plus_reads[row-1,'Chromosome']){
      next
    }
    if (plus_reads[row,'Peak_Start'] <= plus_reads[row-1,'Peak_End']){
      plus_reads[row,'Peak_Start'] <- 
        plus_reads[row-1,'Peak_End']+1 
    }
  }
  for (row in 1:nrow(minus_reads)){
    if (row==nrow(minus_reads)){
      next
    }
    if (minus_reads[row, 'Chromosome'] != minus_reads[row+1,'Chromosome']){
      next
      
    }
    if (minus_reads[row,'Peak_End'] >= minus_reads[row+1,'Peak_Start']){
      minus_reads[row,'Peak_End'] <- minus_reads[row+1,'Peak_Start']-1 
    }
  }
  new_frame <- rbind(plus_reads, minus_reads)
  return(new_frame)
}
head(input_gff)



find_lost_regions <- function (input_gff){
  
  names <- gsub(".*Name= *(.*?) *;.*", "\\1",  input_gff[,9])
  new_gff <- data.frame (input_gff, names)
  only_genes <- new_gff[substring(new_gff[, "names"],1,3) != "id=",]
  list_of_gene_dfs <- split(only_genes, f = list(only_genes$names, only_genes$Chromosome))
  logi <- sapply (list_of_gene_dfs, function(x) nrow(x) >1)
  multiple_peaks <-  list_of_gene_dfs [logi]
  
  lost_regions <- data.frame(multiple_peaks[[1]][1,])
  for (group in  multiple_peaks){

    if (group[1,"Orientation"] == "+"){
      lost_region_start <- group[1,"Peak_End"] 
      lost_region_end <- group[nrow(group),"Peak_End"] 
      if(lost_region_start > lost_region_end){
        print(paste0("skipped",group[1,] ))
        next
      }
      lost_region_row <- group[1,] 
      lost_region_row$Peak_Start <-  lost_region_start 
      lost_region_row$Peak_End <-  lost_region_end
      lost_regions <- rbind(lost_regions, lost_region_row)
      
    }
    else{
      lost_region_start <- group[1,"Peak_Start"] 
      lost_region_end <- group[nrow(group),"Peak_Start"] 
      if(lost_region_start > lost_region_end){
        print(paste0("skipped",group[1,] ))
        next
      }
      lost_region_row <- group[1,] 
      lost_region_row$Peak_Start <-  lost_region_start 
      lost_region_row$Peak_End <-  lost_region_end
      lost_regions <- rbind(lost_regions, lost_region_row)
      
    }
  } 
  return(lost_regions[-1,-10])
}
  
  
write gff <- write.table(x = out, file = "MBCRR_lost_regions.gff", append = F, 
                         quote = F, row.names = F, col.names = F, sep = "\t")
  
  
  
  