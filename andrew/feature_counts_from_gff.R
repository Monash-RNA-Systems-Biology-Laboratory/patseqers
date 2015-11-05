library(varistran)
library(Rsubread)
library(ggplot2)
library(reshape2)
library(stringr)

apa_db_gff <-read.delim("hg19.apadb_v2_final.gff", header=F, comment.char="#",stringsAsFactors=F)

tt_gff<- read.delim("MBCRR-peaks.gff", header=FALSE,
                    comment.char="",stringsAsFactors=F)

annotatum <- data.frame(GeneID = apa_db_gff[,3],
                        Chr = apa_db_gff[,1],
                        Start =apa_db_gff[,4],
                        End = apa_db_gff[,5],
                        Strand = apa_db_gff[,7],
                        stringsAsFactors=FALSE)
feat_counts <- featureCounts(bam_file,annot.ext=annotatum, strandSpecific = 1, nthreads = 8,files = 
                               list.files(pattern = "*.bam$"))

pos_3p <- apply ( 
  feat_counts$annotation, 1, function(row) {
    if (row[5] == "-"){
      return(row[3])
    }
    else if (row[5] == "+"){
      return(row[4])
    }
    else{
      return(NA)
    }
  } 
)
new_anno_df <- data.frame(feat_counts$annotation, pos_3p)
vst_counts<- varistran::vst(feat_counts$counts)
vst_counts<- data.frame(vst_counts, stringsAsFactors = F)
annotated_df <- data.frame (vst_counts, new_anno_df, stringsAsFactors = F)
melted <- melt (annotated_df[,-18])
melted[, 3] <- as.numeric(melted[, 3])
melted[, 4] <- as.numeric(melted[, 4])
names  <- data.frame(t(sapply(melted[,"GeneID"], function(y) strsplit(y,split="\\.")[[1]])))
start_df <- cbind(names, melted[,2:ncol(melted)])
colnames(start_df) <- c("GeneID","Peak", "Chr", "Start", "End", "Strand", "pos_3p", "Sample", "Vst_count")
left_frame <- unique(start_df[, c("GeneID", "Sample")])
split_by_peaks <- split(start_df, start_df$Peak)
counter <- 1
for (peak in split_by_peaks){
  peak$pos_3p <- as.numeric(peak$pos_3p)
  
  new_frame <- peak[, c("GeneID", "Sample", "Vst_count", "pos_3p", "Strand")]
  colnames(new_frame) <-  c("GeneID", "Sample", paste0("peak", counter), paste0("pos_3_", "peak", counter), "Strand")
  left_frame <- merge(left_frame, new_frame, by = c("GeneID", "Sample"), all =T)
  counter <- counter+1
}
# fullish_frame <- left_frame
# fullish_frame$dist <- apply(fullish_frame[,grep("pos_3*", colnames(fullish_frame))], 1, function(row){min(row, na.rm = T)-max(row, na.rm = T)})
trial <- left_frame


melted_peak<- melt(trial, id=c("GeneID","Sample", "Strand.x"), 
                   measure.vars=grep("^peak*", colnames(trial)))

melted_pos <- melt(trial, id=c("GeneID","Sample", "Strand.x"), 
                   measure.vars=grep("^pos*", colnames(trial)))

finished_frame <- data.frame(melted_peak, melted_pos[,3:4])
colnames(finished_frame) <- c("GeneID", "Sample", "Peak", "Vst_count", "pos_3_anno", "pos_3_pos")
#ggplot(data = finished_frame, aes(x= Vst_count, y = pos_3_pos))+geom_point()+facet_wrap(~Sample)
left_frame[is.na(left_frame)] <- 0
finished_frame[is.na(finished_frame)] <- 0


