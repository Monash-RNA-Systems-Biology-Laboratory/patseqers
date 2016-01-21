find_file <- function(file_path) {
    df <- list.files(paste0(file_path, "/expression/genewise/"), pattern = 'counts.csv', recursive = T)
    df <- read.grouped.table(paste0(file_path, "/expression/genewise/", df))
}