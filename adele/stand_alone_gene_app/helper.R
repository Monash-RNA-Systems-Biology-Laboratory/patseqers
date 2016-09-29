all_folders <- as.list(list.dirs(full.names = F, recursive = F))
experiment_list <- all_folders[!grepl("www", all_folders)]

# column_names <- as.list(c("N2.rep1", "N2.rep2", "N2.rep3", "Gld2.rep1", "Gld2.rep2", 
#                   "Gld2.rep3", "Cpb3.rep1", "Cpb3.rep2", "Cpb3.rep3", "Gld2.pcb19.rep1", 
#                   "Gld2.pcb19.rep2", "Gld2.PARN1.rep1", "Gld2.PARN1.rep2", 
#                   "Gld2.CAF1.rep1", "Gld2.CAF1.rep2", "Gld2.PANL2.rep1", 
#                   "Gld2.PANL2.rep2", "Gld2.CCF1.rep1", "Gld2.CCF1.rep2", "Gld2.CCR4.rep1", 
#                   "Gld2.CCR4.rep2", "X1.2.cell.egg.rep1", "X1.2.cell.egg.rep2"))

specific_samples <- as.list(c("N2.rep1", "N2.rep2", "N2.rep3", "Gld2.rep1", "Gld2.rep2", 
                          "Gld2.rep3", "Cpb3.rep1", "Cpb3.rep2", "Cpb3.rep3", "Gld2.pcb19.rep1", 
                          "Gld2.pcb19.rep2", "Gld2.CCF1.rep1", "Gld2.CCF1.rep2", "X1.2.cell.egg.rep1", "X1.2.cell.egg.rep2"))


find_files <- function(file_path) {    
    count_df <- list.files(paste (file_path), pattern = 'genewise-count.csv', recursive = T)
    #columnNames <- read.csv(paste0(file_path, "/", count_df), header = F, row.names = 1, nrow = 1)
    count_df <- read.csv(paste0(file_path, "/", count_df), row.names = 1)
    
    tail_df <- list.files(paste (file_path), pattern = 'genewise-tail.csv', recursive = T)
    tail_df <- read.csv(paste0(file_path, "/", tail_df), row.names = 1)
    
    tail_count_df <- list.files(paste (file_path), pattern = 'genewise-tail-count.csv', recursive = T)
    tail_count_df <- read.csv(paste0(file_path, "/", tail_count_df), row.names = 1)
    
    info_df <- list.files(paste (file_path), pattern = 'genewise-info.csv', recursive = T)
    info_df <- read.csv(paste0(file_path, "/", info_df), row.names = 1)
    
    all_files <- list(count_df, tail_df, tail_count_df, info_df)
    return(all_files)
}

### Normalise data

tmmnorm <- function(data) {
    b <- DGEList(data)
    b <- calcNormFactors(b)
    b <- cpm(b, log=TRUE, prior.count = 0.5)
}

#The first arguement "attr" chooses the gene nonmenclature to be returned from the query. 
#This is dependentant on the organism and dataset. The Ce.elegans data uses refseq identifies whereas the
#yeast test data uses ensembl gene ids. 'pick.mart' is the organism mart to used in the query. Term is the
#input GOterm the query will use.

fetch_goterm <- function(attr, term, pick_mart) {
    get.GO <- getBM(attributes=c(attr), filters = "go_id",
                    values = term, mart= pick_mart)
    return(get.GO)
}

# fetch_names <- function(pick_mart, term) {
#     get.names <- getBM(attributes=c("refseq_mrna"), filters = "go_id",
#                     values = term, mart= pick_mart)
#     return(get.GO)
# }



#Create GO-Term
#ensem <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "celegans_gene_ensembl", host = "www.ensembl.org")
#attr <- "refseq_mrna"

#all_GOslim <- getBM(attributes=c("goslim_goa_accession"), mart= ensem)
#write.csv(all_GOslim, "GOslim_accessions.txt")


goslim.imp <- read.csv("GOslim_accessions.txt")
goslim.imp$X <- NULL
colnames(goslim.imp) <- "GOSlim"
goslim_list <- as.list(as.vector(goslim.imp$GOSlim))

