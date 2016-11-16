## Define the list for the 
specific_samples <- as.list(c("N2.rep1", "N2.rep2", "N2.rep3", "Gld2.rep1", "Gld2.rep2", 
                          "Gld2.rep3", "Cpb3.rep1", "Cpb3.rep2", "Cpb3.rep3", "Gld2.pcb19.rep1", 
                          "Gld2.pcb19.rep2", "Gld2.CCF1.rep1", "Gld2.CCF1.rep2", "X1.2.cell.egg.rep1", "X1.2.cell.egg.rep2"))

## Pull out the files of interest from the gld2 pipeline 
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


### RDS objects to make the app run faster: 

## These steps can be run from the console seperately. 
## To define the RDS object for the GO Slim Accession list that'll appear in a select input:

## Create an C elegans ensembl mart object
#marty.mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "celegans_gene_ensembl", host = "www.ensembl.org")

## Return both GO Slim IDs and their description
#goslims_both <- getBM(attributes=c("goslim_goa_accession", "goslim_goa_description"), mart= marty.mart)

## Create a named list and loop through the dataframe - names are the GO Slim ID and name, with the GO Slim ID alone as the item to the list 
# slim_list <- list()
# for (i in 1:nrow(goslims_both)) {
#   name_list_item <- paste0(goslims_both$goslim_goa_description[i], " (" , goslims_both$goslim_goa_accession[i], ")")
#   slim_list[name_list_item] <- goslims_both$goslim_goa_accession[i]
# }

## Save the list as an RDS
#save(slim_list, file = "www/data/preload_goslim_list.rds")

## These steps need to be run from server only once, they'll generate the starting objects to be used in server. 

# Call count, tail, tail count and info files from the gld2 pipeline directory 
#file_call <- find_files("gld2combined")

# Save the list of dataframes as an RDS object
#save(file_call, file = "www/data/preload_files_called.rds")

# Create a C.elegans ensembl mart object. 
#selected.mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "celegans_gene_ensembl", host = "www.ensembl.org")

# Save the mart object as an RDS object.
#save(selected.mart, file = "www/data/preload_mart.rds")