library(edgeR)
library(biomaRt)

experiment_list <- as.list(list.dirs(full.names = F, recursive = F))

find_files <- function(file_path) {    
    count_df <- list.files(paste (file_path), pattern = 'genewise-count.csv')
    count_df <- read.csv(paste0(file_path, "/", count_df), row.names = 1)
    
    tail_df <- list.files(paste (file_path), pattern = 'genewise-tail.csv')
    tail_df <- read.csv(paste0(file_path, "/", tail_df), row.names = 1)
    
    tail_count_df <- list.files(paste (file_path), pattern = 'genewise-tail-count.csv')
    tail_count_df <- read.csv(paste0(file_path, "/", tail_count_df), row.names = 1)
    
    info_df <- list.files(paste (file_path), pattern = 'genewise-info.csv')
    info_df <- read.csv(paste0(file_path, "/", info_df), row.names = 1)
    
    all_files <- list(count_df, tail_df, tail_count_df, info_df)
    return(all_files)
}

tmmnorm <- function(data) {
    b <- DGEList(data)
    b <- calcNormFactors(b)
    b <- cpm(b, log=TRUE, prior.count = 0.5)
}

#The first arguement "attr" chooses the gene nonmenclature to be returned from the query. 
#This is dependentant on the organism and dataset. The Ce.elegans data uses refseq identifies whereas the
#yeast test data uses ensembl gene ids. 'pick.mart' is the organism mart to used in the query. Term is the
#input GOterm the query will use.

fetch_goterm <- function(attr, pick_mart, term) {
    get.GO <- getBM(attributes=c(attr), filters = "go_id",
                    values = term, mart= pick_mart)
    return(get.GO)
}

# fetch_names <- function(pick_mart, term) {
#     get.names <- getBM(attributes=c("refseq_mrna"), filters = "go_id",
#                     values = term, mart= pick_mart)
#     return(get.GO)
# }

find_peaks <- function(file_path) {    
    peaks_df <- list.files(paste (file_path), pattern = 'individual-pairs.csv')
    peaks_df <- read.grouped.table(paste0(file_path, "/", count_df))
}