
consistent_up_down <- function(input_csv){
  dg_csv <- read.csv (input_csv)
  no_hla <- dg_csv[!grepl(pattern = "*HLA-*", x =dg_csv$gene),]
  
  all_up <-no_hla[no_hla$LM2 < no_hla$HM & no_hla$LM2 >0, ]
  all_down <- no_hla[no_hla$LM2 > no_hla$HM & no_hla$LM2 < 0, ]
  ordered_up <- all_up [with (all_up, order (HM, LM2, decreasing = T)),]
  ordered_down <- all_down  [with (all_down , order (HM, LM2 )),]
  
}


