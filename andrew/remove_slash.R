slash_fixer <- function (slashed_name_file) {
  input <- read.delim(slashed_name_file, header = T,stringsAsFactors=FALSE )
  input <- input [,1:2]
  colnames(input)<- c("gene", "FDR")
  new_frame <-data.frame(       gene=character(), 
                                FDR=character(), 
                                stringsAsFactors=FALSE) 
  for (i in 1:nrow(input)){
    split <- strsplit(input[i,1] , "/")
    if (length(split[[1]])> 1 ){
      for (name in split[[1]][2:length(split)]){
          new_name <-  data.frame(name,input[i,2])
          colnames(new_name)<- c("gene", "FDR")
          new_frame <- rbind(new_frame,new_name) 
       
      }    
      
    }
    else{
      new_frame <- rbind(new_frame,input[i,]) 
    }
    
  }
  write.table(x=new_frame, file = paste0("unslashed", slashed_name_file), sep = "\t",quote = F, row.names = F, col.names = F)
  return(new_frame)
}