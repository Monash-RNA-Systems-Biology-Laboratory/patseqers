#' @export
make_o <-
function(mat, enable=TRUE, fast=FALSE) {
    # Note: paths are given ordered by order
    mat[is.na(mat)] = 0.001
    if (nrow(mat) < 3 || !enable) {
        list(
            dendrogram = NULL,
            order = seq_len(nrow(mat)),
            paths = rep('',nrow(mat)))
        
    } else if (fast) {
        # Too many for seriation to deal with efficiently
        perm <- seriation::seriate(mat, method='PCA')[[1]]
        
        list(dendrogram = NULL,
             order = seriation::get_order(perm),
             paths = rep('',nrow(mat)))
        
    } else {
        dist_mat <- dist(mat)
        control <- list(hclust = hclust(dist_mat))
        dend_mat <- as.dendrogram(
            seriation::seriate(dist_mat,
                               method = 'OLO',
                               control = control)[[1]])
        
        list(dendrogram = dend_mat,
             order = order.dendrogram(dend_mat),
             paths = dendrogram_paths(dend_mat))
    }
}
