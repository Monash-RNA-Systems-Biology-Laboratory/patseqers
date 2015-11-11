# Define some gene sets. 
library(ggplot2)
library(reshape2)
consistent_up_down <- function(input_csv){
  dg_csv <- read.csv (input_csv)
  only_signif_fdr <- dg_csv[dg_csv$FDR <0.05, ]
  no_hla <- only_signif_fdr[!grepl(pattern = "*HLA-*", x =only_signif_fdr$gene),]
  
  all_up <-no_hla[no_hla$LM2 < no_hla$HM & no_hla$LM2 >0, ]
  all_down <- no_hla[no_hla$LM2 > no_hla$HM & no_hla$LM2 < 0, ]
  ordered_up <- all_up [with (all_up, order (HM, LM2, decreasing = T)),]
  ordered_down <- all_down  [with (all_down , order (HM, LM2 )),]
  
}

lung_homing <- function (input_csv){
  dg_csv <- read.csv (input_csv)
  only_signif_fdr <- dg_csv[dg_csv$FDR <0.05, ]
  no_hla <- only_signif_fdr[!grepl(pattern = "*HLA-*", x =only_signif_fdr$gene),]
  
  lung_genes <-no_hla[no_hla$LM2 > no_hla$HM & 
                        no_hla$LM2 >0 & no_hla$HM < no_hla$LM2 & no_hla$LM2 > 3, ]
  ordered_up <-  lung_genes [with ( lung_genes, order (LM2,HM, decreasing = T)),]
  return(ordered_up)
}



write.table(x = ordered_up, file = "MBCR_all_lines_consistent_up.csv", append = F, 
                         quote = F, row.names = F, col.names = F, sep = ",")

write.table(x = ordered_down, file = "MBCR_all_lines_consistent_down.csv", append = F, 
            quote = F, row.names = F, col.names = F, sep = ",")
ggplot(data = ordered_down, aes( x= ,y = HM))



human_3_p_cpsf_complex <- list ("CPSF1", "CPSF2", "CPSF3", "CPSF4", "hFip1")
