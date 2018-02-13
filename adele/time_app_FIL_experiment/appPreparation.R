### This contains the pre-processing work to create the data used by the app

# Read in the grouped.table counts file (csv that contains different outputs from the tail tools pipeline). See this link for more: https://github.com/Victorian-Bioinformatics-Consortium/tail-tools/blob/master/doc/statistics.md
# The individual csvs could be used but here this just allows both the annotation data frame and the counts data frame to be loaded in one step
# data_imp <- nesoni::read.grouped.table("counts.csv")

# Save to RDS
# save(data_imp, file = "dat.RDS")

### Load the count.csv that was produced by the tail tools pipeline. This is a nesoni grouped table
load("dat.RDS")
    
#Pull out the count matrix and the annotations for the count matrix
count_mat_preprocessing <- data_imp$Count
annot_df <- data_imp$Annotation

# Pull out the C.albicans systematic name and the standard name
annot_systematic_name <- rownames(annot_df)
annot_standard_name <- annot_df$gene

# Create an empty vector then loop through the systematic names and paste them together with the standard names. This creates a vector containing both within one string, separated with a '/'. This allows either annotation to be searched from the drop down select menu.
gene_names <- c()
for(i in 1:length(annot_systematic_name)) {
  gene_names[i] <- paste0(annot_standard_name[i], "/", annot_systematic_name[i])
}

# Assign the new combined names as the row names to the count data set
rownames(count_mat_preprocessing) <- gene_names

# Rename the columns to get the file name conventions between the two replicate sets identical
colnames(count_mat_preprocessing) <- gsub("^FIL_", "R2_FIL_", colnames(count_mat_preprocessing))
colnames(count_mat_preprocessing) <- gsub("^INHIB_", "R2_INHIB_", colnames(count_mat_preprocessing))
colnames(count_mat_preprocessing) <- gsub("^Minus", "R1_FIL_", colnames(count_mat_preprocessing))
colnames(count_mat_preprocessing) <- gsub("^Plus", "R1_INHIB_", colnames(count_mat_preprocessing))
colnames(count_mat_preprocessing) <- gsub("^YPD_0", "B1_YPD_R", colnames(count_mat_preprocessing))
colnames(count_mat_preprocessing) <- gsub("^YPD_R", "B2_YPD_R", colnames(count_mat_preprocessing))

# Save to an RDS object. This will be loaded into the app
#save(count_mat_preprocessing, file = "proprocessed_counts.RDS")
