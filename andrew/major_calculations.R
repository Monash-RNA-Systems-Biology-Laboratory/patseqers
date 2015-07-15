#File paths in place for when Michael directs me. 
gff_file_path <- "/home/bigpatto2/Bioinformatics/SCP_2.0/SHINY_yeast_test/"
bam_file_path <- "/home/bigpatto2/Bioinformatics/SCP_2.0/SHINY_yeast_test/"

#List all datasets inside a directories with "SHINY" in the name. 

list_datasets <- function(){
  return(c("Yeast", "Cancer"))
}

# Returns a list of bam files from the nominated directory
find_bam_files <- function(file_path) {
  bam_files <- list.files(paste (file_path), pattern = '*.bam$')
  return(bam_files)
}
# Returns a list of gff files from the nominated directory
find_gff_files <- function(file_path) {
  gff_files <- list.files(paste(file_path), pattern = '*.gff$')
  return(gff_files)
}

make_plot <- function(processed_frame, ranges,names, leg,group){
  if (group == T){
  samples <- split(processed_frame, processed_frame$group, drop =T)
  }
  else {
    samples <- split(processed_frame, processed_frame$sample, drop =T)    
  }  
  par(bty="l")
  par(mar=c(5.1,4.1,4.1,8.1), xpd =T)
    
  dummy_ecdf <- ecdf(1:10)
  curve((-1*dummy_ecdf(x)*100)+100, from=ranges[1], to=ranges[2], 
        col="white", xlim=ranges, main= paste(names),
        axes=F, xlab= 'Poly (A)-Tail Length', ylab = 'Percent Population (%)', ylim =c(0,100))
  axis(1, pos=0, tick = 25)
  axis(2, pos= 0, at= c(0,25,50,75,100), tick = 25)
  
  
  count <- 1
  
  
  for (df in samples){
   split_peak <- split(df,df$gene_or_peak_name, drop =T)    
   for(gene_or_peak in split_peak){  
      colours <- rainbow(length(samples)*length(split_peak))
      ecdf_a <- ecdf(gene_or_peak[,"number_of_as"])
      curve((-1*ecdf_a(x)*100)+100, from=ranges[1], to=ranges[2], 
                     col=colours[count], xlim=ranges, main= paste(names),
                     add=T)
      count <- count +1     
        
  }
  
  }
  leg_names <- list()
  for (name in names(samples)){
    leg_names <-c(leg_names, paste(name, names(split_peak)))
    
  }
  
  if (leg ==T){ 
    legend(ranges[2]-40,95 + (length(samples)*5), 
           legend = leg_names, fill = colours, bty ="n")
  }

}

selected_data<- function (data){
  return (paste(data))
}
# Outpus the rows matching the input gene or peak name   
filter_gff_for_rows<- function (gff,names){
    split_names <- strsplit(names, split = " ")
    empty <- data.frame()
    
    for (name in split_names[[1]]){
      index1 <- with(gff, grepl 
                     (ignore.case = T,paste('=',name,'[;/$]{1}',sep=""), gff[,'Information']))
      # Would be nice to find some better regex to get rid of this if statement. 
      # Maybe do this with a GFF parser
      
      if (sum(index1)==0){
        index1 <- with(gff, grepl 
                       (ignore.case = T,paste('=',name,'$',sep=""), gff[,'Information']))
      }
      
      output <-gff[index1, ] 
      empty <- rbind(empty, output)
    }
    return(empty)
}

# This function gets the poly (A) counts for all given gff rows
get_a_counts <- function(bam_file_path,gff_rows, bam_files,names, groups){
  reads_report <- data.frame() 
  split_names <- strsplit(names, " ")
  
  for (gff_row in 1:nrow(gff_rows)){
    counts_frame <- get_a_counts_gff_row(bam_file_path, gff_rows[gff_row,], 
                                         bam_files, groups)
    counts_frame$gene_or_peak_name <- split_names[[1]][gff_row]
    reads_report <-rbind(reads_report,counts_frame)      
  }
  return(reads_report)
}


get_a_counts_gff_row <- function(bam_file_path,peak, bam_files, groups){
  if (peak[,"Orientation"]== "-"){
    ori <- TRUE    
  }
  else{
    ori <- FALSE
  }
  bam_frame <- data.frame()
  count <- 1 
  for (bam_file in bam_files){
    full_file_path <-paste(bam_file_path,bam_file, sep ="")
        
    param <- ScanBamParam(what=c('qname','pos','qwidth','strand'),
                          tag=c('AN','AD'),flag=scanBamFlag(isMinusStrand=ori) , 
                          which=GRanges(peak [,'Chromosome'],IRanges(
                          peak[,'Peak_Start'], peak[,'Peak_End'] )))
    #Grabs reads overlapping the range specified by the gff row
    result <- scanBam (full_file_path , param = param, isMinusStrand = ori)
    # A check to make sure the adapter bases column is present. 
    #If not, I make a fake one of NAs.
    if (length(result [[1]][[5]][[2]])!= length(result [[1]][[5]][[1]])){
      result [[1]][[5]][[2]] <- rep(NA, length(result [[1]][[5]][[1]]))      
    }
    
    single_bam_frame <-  data.frame(result) 
    colnames(single_bam_frame)<- c("qname", "strand", "pos", 
                                   "width", "number_of_as", "number_of_ad_bases")
    #If the read is on the forward strand, add width to pos to obtain 3' end. 
    if (ori == FALSE ){
      single_bam_frame$pos <- single_bam_frame$pos+ single_bam_frame$width
    }
    single_bam_frame <- single_bam_frame[single_bam_frame$pos >= 
                                           peak[,'Peak_Start']& 
                                           single_bam_frame$pos <= peak[,'Peak_End'] ,]
    single_bam_frame$sample <- paste(bam_file)
    single_bam_frame$group<- paste("group", groups[count])
    bam_frame <- rbind(bam_frame,single_bam_frame)
    count <- count +1
    
  }  
  return(bam_frame)
}

#Strips whitespace out of names
trim <- function (x) gsub("^\\s+|\\s+$", "", x)


names_string <- function(s_frame, groups){
  to_print <- character()
  for (frame in s_frame){
    split_peaks <- split(frame ,frame$gene_or_peak_name, drop =T)
    for (peak_frame in split_peaks){
        
      if (groups == T){
        str <- paste("The poly (A) read count for ", peak_frame$group[1]," ",
                     peak_frame$gene_or_peak_name[1], " is: ",nrow(peak_frame),".", "\n", sep ="")
      }
      else{
        str <- paste("The poly (A) read count for ", peak_frame$sample[1]," ",
                     peak_frame$gene_or_peak_name[1], " is: ",nrow(peak_frame),".", "\n", sep ="")
    }
    to_print <- c(to_print, str)
    }
  }
  return(to_print)
}
modify_gff_inplace <- function (gff_file) {
  
  start_gff_file <- read.delim(gff_file, header=FALSE,
                               comment.char="",stringsAsFactors=F)
  start_gff_file<- start_gff_file[-1,]
  colnames(start_gff_file)<- c('Chromosome', 'Generated_By', 'Feature_Type', 
                               'Peak_Start','Peak_End','-',
                               'Orientation', '--','Information')
  new_frame <- start_gff_file[
    with(start_gff_file,order(
      Chromosome,Orientation,Peak_Start)
    ),
    ]
  for (row in 1:nrow(new_frame)){
    
    if (new_frame[row,'Orientation'] == '+'){
      new_frame[row,c('Peak_End', 'Peak_Start')] <- 
        new_frame[row,c('Peak_End', 'Peak_Start')]+10
      if (new_frame[row,'Peak_Start'] <=
            new_frame[row-1,'Peak_End']&
            new_frame[row,'Chromosome'] ==
            new_frame[row-1,'Chromosome']){
        new_frame[row,'Peak_Start'] <- 
          new_frame[row-1,'Peak_End']+1          
      }
    }
    else{
      if (new_frame[row,'Peak_End'] >=
            new_frame[row+1,'Peak_Start']&
            new_frame[row, 'Chromosome'] ==
            new_frame[row+1,'Chromosome']){
        new_frame[row,'Peak_End'] <- 
          new_frame[row+1,'Peak_Start']-1          
      }
      new_frame[row,c('Peak_End', 'Peak_Start')] <- 
        new_frame[row,c('Peak_End', 'Peak_Start')]-10
    }
  }
  return(new_frame)
}