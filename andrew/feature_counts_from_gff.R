library(varistran)
library(Rsubread)
library(ggplot2)

apa_db_gff <-read.delim("hg19.apadb_v2_final.gff", header=F, comment.char="#",stringsAsFactors=F)

tt_gff<- read.delim("MBCRR-peaks.gff", header=FALSE,
                    comment.char="",stringsAsFactors=F)

annotatum <- data.frame(GeneID = apa_db_gff[,3],
                        Chr = apa_db_gff[,1],
                        Start =apa_db_gff[,4],
                        End = apa_db_gff[,5],
                        Strand = apa_db_gff[,7],
                        stringsAsFactors=FALSE)

feat_counts_df <- data.frame()
for (bam_file in list.files(pattern = "*.bam$")){
  feat_counts <- featureCounts(bam_file,annot.ext=annotatum)
  df <- data.frame(feat_counts$counts)
  if (length(feat_counts_df) == 0){
    feat_counts_df <- df
  }
  else{
    feat_counts_df <- data.frame(feat_counts_df, df)    
  }   
}
feat_counts_df <- data.frame(varistran::vst(feat_counts_df))

# ggplot(data =feat_counts_df, aes(x = LNA_rep_1.bam, y= MDA_MB_231_sorted.bam) )+
#   geom_point()



