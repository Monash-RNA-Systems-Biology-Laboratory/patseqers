
# Returns a list of bam files from the nominated directory
find_bam_files <- function(file_path) {
  if(file.exists(paste0(file_path,'/','plotter-config.json'))){
    json <- fromJSON(paste0(file_path,'/', "plotter-config.json"))
    bam_files <- json$samples
  }
  else{
    bam_files <- list.files(paste (file_path), pattern = '*.bam$')
  }
  return(bam_files)
}

# Returns a list of gff files from the nominated directory
find_gff_files <- function(file_path) {
  if(file.exists(paste0(file_path,'/','plotter-config.json'))){
    json <- fromJSON(paste0(file_path,'/', "plotter-config.json"))
    gff_files <- json$peaks    
  }
  else{
    gff_files <- list.files(paste(file_path), pattern = '*.gff$')    
  }
  return(gff_files)
}
# So this function got out of hand... but I basically makes whatever plot
# is specified by the user. 
make_plot <- function(processed_frame, ranges,names, leg,group, alt_plot, order_alt, alt_cumu_dis, poly_a_pileup,show_poly_a =F){
  if (length(processed_frame) == 0){
    return("Error, no reads for this gene or eak in this sample")
  }
  
  if(order_alt==T){

    new_frame <- processed_frame[
      with(processed_frame,order(
        -width, -number_of_as)
      ),
      ]
    ylab <- "Sorted Read Number"
  }
  else{
    new_frame <- processed_frame
    ylab <- "Read Number"
  }
  if (group == T){
    samples <- split(new_frame, new_frame$group, drop =T)
  }
  else {
    samples <- split(new_frame, new_frame$sample, drop =T)    
  }  
 
  par(bty="l", ps = 10, mar=c(5.1,4.1,4.1,8.1), xpd =T)
  if(alt_plot == T){
    
    if (poly_a_pileup == T ){
      if (length(samples) == 1){
        par(mfrow= c(1,1))
      }
      else if ((length(samples)/2)%%1 == 0){
        par(mfrow= c(as.integer(length(samples)/2),2))
      }
      else{
        par(mfrow= c(as.integer(length(samples)/2)+1,2))        
      }
      for (sample in samples) {
        points <- data.frame(sample$width, sample$number_of_as)
        ymax <- nrow(points)  
        
        count <- 1:ymax
        
        plot(NA,xlim=ranges, ylim = c(0, ymax), xlab= "Number of Bases", ylab = ylab, 
             main= paste(sample[1,'sample']))
        for (i in 1:ymax){
          segments(x0= 0, y0= i,x1= points[i,1], col="purple")
          segments(x0= points[i,1], y0= i,x1= points[i,1] +points[i,2] , col="pink")
          
        }
        
      }
      return()
    }
    seq_to_plot <- get_plot_sequence(new_frame)
    ymax <- 0
    for (sample in samples){
      title <- sample[1, 'gene_or_peak_name']
      if (nrow (sample) > ymax){
        ymax <- nrow(sample)
      }
      
    }
    
    if (alt_cumu_dis ==T) {
      dummy_ecdf <- ecdf(1:10)
      curve((-1*dummy_ecdf(x)*100)+100, from=ranges[1], to=ranges[2], 
            col="white", xlim=ranges, main= paste(names),
            axes=F, xlab= "Number of Bases", ylab = 'Percent Population (%)', ylim =c(0,100))
      axis(1, pos=0, tick = 25)
      axis(2, pos= 0, at= c(0,25,50,75,100), tick = 25) 
      count <- 1  
      for (df in samples){
        split_peak <- split(df,df$gene_or_peak_name, drop =T)    
        for(gene_or_peak in split_peak){  
          colours <- rainbow(length(samples)*length(split_peak))
          
          ecdf_a <- ecdf(gene_or_peak[,"width"])
          curve((-1*ecdf_a(x)*100)+100, from=ranges[1], to=ranges[2], 
                col=colours[count], xlim=ranges, main= paste(names),
                add=T)
          count <- count +1     
          
        }
        # This loop makes a list for the legend. 
        leg_names <- list()
        for (name in names(samples)){
          leg_names <- c(leg_names, paste(name, names(split_peak)))
          
        }
        if (leg ==T){ 
          x_offset <-  length(strsplit(paste(leg_names), "")[[1]])
          legend(ranges[2]-30-(x_offset)*2,110 +(length(samples)*-0.8), 
                 legend = leg_names, fill = colours, bty ="n")
        }
        
      }
      
    }
    else{
      plot(NA,xlim=ranges, ylim = c(0, ymax), xlab= "Number of Bases", ylab = ylab, 
           main= title)
      
      
      col_count <- 1
      for (df in samples){
        split_peak <- split(df,df$gene_or_peak_name, drop =T)    
        colours <- rainbow(length(samples)*length(split_peak))
        for (sample in split_peak) {
          
          points <- data.frame(sample$width, sample$number_of_as) 
          count <- 1:nrow(points)
          
          lines(points[,1],count,col = colours[col_count])
          
          if (show_poly_a ==T){
            lines(points[,1]+ points[,2], count, col = colours[col_count])  
          }
          col_count <- col_count +1
        }
        
        leg_names <- list()
        for (name in names(samples)){
          leg_names <- c(leg_names, paste(name, names(split_peak)))
          
        }
        if (leg ==T){ 
          legend(y = ymax+38,x = ranges[2]-55, 
                 legend = leg_names, fill = colours, bty ="n")
        }
        
      }
      
    }
    #single_bases <- strsplit(seq_to_plot, "")
    #axis(1, pos=-7,at = 1:nchar(seq_to_plot), labels = single_bases[[1]], lwd =0, cex.axis = 0.3) 
    
  }
  
  else{    
    
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
    # This loop makes a list for the legend. 
    leg_names <- list()
    for (name in names(samples)){
      leg_names <- c(leg_names, paste(name, names(split_peak)))
      
    }
    if (leg ==T){ 
      x_offset <-  length(strsplit(paste(leg_names), "")[[1]])
      legend(ranges[2]-30-(x_offset)*2,110 +(length(samples)*-0.8), 
             legend = leg_names, fill = colours, bty ="n")
    }
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
                   (ignore.case = T,paste('[=/]{1}',name,'[;/,]',sep=""), gff[,'Information']))
    # Would be nice to find some better regex to get rid of this if statement. 
    # Maybe do this with a GFF parser
    
    index2 <- with(gff, grepl 
                   (ignore.case = T,paste('=',name,'$',sep=""), gff[,'Information']))
    
    output <-gff[index1 | index2, ] 
    if (nrow (output) == 0){
      stop('There are no reads this gene/peak in your selected samples')
    }
    output$input_gene_or_peak <- name
    empty <- rbind(empty, output)
  }

  return(empty)
}

# This function gets the poly (A) counts for all given gff rows
get_a_counts <- function(bam_file_path,gff_rows, bam_files, groups, names_from_json){
  reads_report <- data.frame() 
  for (gff_row in 1:nrow(gff_rows)){
    print("here?")
    counts_frame <- get_a_counts_gff_row(bam_file_path, gff_rows[gff_row,], 
                                         bam_files, groups, names_from_json)
    print("shouldn't get here")
    if (nrow(counts_frame) == 0){
      next
    }
    counts_frame$gene_or_peak_name <- gff_rows[gff_row, 'input_gene_or_peak']

    reads_report <-rbind(reads_report,counts_frame)   

  }
  
  return(reads_report)
}

# Parses the BAM files for eah GFF file entry that we are given
get_a_counts_gff_row <- function(bam_file_path,peak, bam_files, groups,names_from_json){
  if (peak[,"Orientation"]== "-"){
    ori <- TRUE    
  }
  else{
    ori <- FALSE
  }
  bam_frame <- data.frame()
  count <- 1
  for (bam_file in bam_files){
    if (substring(bam_file,1,1)=="/"){
      full_file_path <-bam_file
    }
    else{
      full_file_path <-paste(bam_file_path,"/", bam_file, sep ="")
    }
    
    param <- ScanBamParam(what=c('qname','pos','qwidth','strand', 'seq'),
                          tag=c('AN','AD'),flag=scanBamFlag(isMinusStrand=ori) , 
                          which=GRanges(peak [,'Chromosome'],IRanges(
                            peak[,'Peak_Start'], peak[,'Peak_End'] )))
    #Grabs reads overlapping the range specified by the gff row
    result <- scanBam (full_file_path , param = param, isMinusStrand = ori)
    # A check to make sure the adapter bases column is present. 
    #If not, I make a fake one of 0s.
    
    if (length(result [[1]][[6]][[1]])!= length(result [[1]][[5]])){
      result [[1]][[6]][[1]] <- rep(0, length(result [[1]][[5]]))    
    }
    if (length(result [[1]][[6]][[2]])!= length(result [[1]][[5]])){
      result [[1]][[6]][[2]] <- rep(0, length(result [[1]][[5]]))      
    }
    result[[1]][["seq"]] <- as.character(result[[1]][["seq"]])
    if (length(result [[1]][[5]]) == 0){
      stop(paste('There are no reads for at least one peak in ', bam_file))
    }
    
    single_bam_frame <-  data.frame(result) 
    
    colnames(single_bam_frame)<- c("qname", "strand", "pos", 
                                   "width", "sequence", "number_of_as", "number_of_ad_bases")
    #If the read is on the forward strand, add width to pos to obtain 3' end. 
    if (ori == FALSE ){
      single_bam_frame$pos <- single_bam_frame$pos+ single_bam_frame$width
    }
    
    single_bam_frame <- single_bam_frame[single_bam_frame$pos >= 
                                           peak[,'Peak_Start']& 
                                           single_bam_frame$pos <= peak[,'Peak_End'] ,]
    print(str(single_bam_frame))
    if (nrow(single_bam_frame) == 0){
      next
    #  single_bam_frame <- data.frame(0,0,0,0,0,0,0,0,0,0)
    }
    if (substring(bam_file,1,1)=="/"){
      single_bam_frame$sample <-  names_from_json$name [names_from_json$bam ==paste(bam_file)]
    }
    else{
      single_bam_frame$sample <- paste(bam_file)
    }
    print("do we get to here?")
    single_bam_frame$group<- paste("group", groups[count])
 
    bam_frame <- rbind(bam_frame,single_bam_frame)
    count <- count +1
    
  }  
  return(bam_frame)
}

#Strips whitespace out of names
trim <- function (x) gsub("^\\s+|\\s+$", "", x)  

# 
names_string <- function(s_frame, groups, all_reads){
  # The data frame is split by samples here
  to_print <- character()
  for (frame in s_frame){
    #The data frame is split into genes here
    split_peaks <- split(frame ,frame$gene_or_peak_name, drop =T)
    for (peak_frame in split_peaks){
      if (all_reads == T){
        tail_reads <- " "        
      }
      else{
        tail_reads <- " with a poly (A)-tail "
      }      
      if (groups == T){
        str <- paste("The number of reads ",tail_reads,"for ", peak_frame$group[1]," ", 
                     peak_frame$gene_or_peak_name[1], " is: ",nrow(peak_frame),".", "\n", sep ="")
      }
      else{
        str <- paste("The number of reads",tail_reads,"for ",  peak_frame$sample[1]," ",
                     peak_frame$gene_or_peak_name[1], " is: ",nrow(peak_frame),".", "\n", sep ="")
      }
      to_print <- c(to_print, str)
    }
  }
  return(to_print)
}
# Handles overlapping peaks in the gff file 
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
# Makes the menas and medians frame shown in info tab of the app
make_means_and_meds_frame <- function (poly_a_counts){
  if (poly_a_counts[1,"group"]== "group NULL"){
    into_samples <- split(poly_a_counts, 
                          list(poly_a_counts$sample, poly_a_counts$gene_or_peak_name))
    mm_frame <- data.frame()
    for (sample in into_samples){
      sample_mean <- mean(sample$number_of_as, na.rm =T)
      sample_median <- median(sample$number_of_as, na.rm =T)
      name <- paste(sample[1, "sample"], sample[1, "gene_or_peak_name"])
      to_bind <- cbind(name,sample_mean, sample_median)
      mm_frame <- rbind(mm_frame,to_bind)
    }
    
  }
  else{
    into_samples <- split(poly_a_counts, poly_a_counts$group)
    mm_frame <- data.frame()
    for (sample in into_samples){
      sample_mean <- mean(sample$number_of_as, na.rm =T)
      sample_median <- median(sample$number_of_as, na.rm =T)
      to_bind <- cbind(sample[1, "group"],sample_mean, sample_median)
      mm_frame <- rbind(mm_frame,to_bind)
    }
  }  
  colnames (mm_frame) <- c("Sample Name", "Mean Poly (A)-Tail Length", "Median Poly (A)-Tail Length")
  return(mm_frame)
}




