\name{pl_hmap_detailed}
\alias{pl_hmap_detailed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{pl_hmap_detailed
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pl_hmap_detailed(matf, gmatf, clusterby = 1, cluster_features = 1, sample_labels = NULL, sample_labels2 = NULL, feature_labels = NULL, gene_labels = NULL, product_labels = NULL, feat_span = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matf}{
%%     ~~Describe \code{matf} here~~
}
  \item{gmatf}{
%%     ~~Describe \code{gmatf} here~~
}
  \item{clusterby}{
%%     ~~Describe \code{clusterby} here~~
}
  \item{cluster_features}{
%%     ~~Describe \code{cluster_features} here~~
}
  \item{sample_labels}{
%%     ~~Describe \code{sample_labels} here~~
}
  \item{sample_labels2}{
%%     ~~Describe \code{sample_labels2} here~~
}
  \item{feature_labels}{
%%     ~~Describe \code{feature_labels} here~~
}
  \item{gene_labels}{
%%     ~~Describe \code{gene_labels} here~~
}
  \item{product_labels}{
%%     ~~Describe \code{product_labels} here~~
}
  \item{feat_span}{
%%     ~~Describe \code{feat_span} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (matf, gmatf, clusterby = 1, cluster_features = 1, sample_labels = NULL, 
    sample_labels2 = NULL, feature_labels = NULL, gene_labels = NULL, 
    product_labels = NULL, feat_span = 1) 
{
    matf <- as.data.frame(matf)
    tvec <- grepl("Tail", colnames(matf))
    y <- as.matrix(matf[, tvec])
    z <- as.matrix(matf[, !tvec])
    if (is.null(sample_labels) && !is.null(colnames(y))) 
        sample_labels <- colnames(y)
    if (is.null(sample_labels)) 
        sample_labels <- rep("", ncol(y))
    sample_labels[is.na(sample_labels)] <- ""
    if (is.null(sample_labels2) && !is.null(colnames(z))) 
        sample_labels2 <- colnames(z)
    if (is.null(sample_labels2)) 
        sample_labels2 <- rep("", ncol(z))
    sample_labels[is.na(sample_labels)] <- ""
    if (is.null(feature_labels) && !is.null(rownames(y))) 
        feature_labels <- rownames(y)
    if (is.null(gene_labels) && !is.null(gmatf$gene)) 
        gene_labels <- gmatf$gene
    if (is.null(product_labels) && !is.null(gmatf$product)) 
        product_labels <- gmatf$product
    if (is.null(feature_labels)) 
        feature_labels <- rep("", nrow(y))
    if (is.null(gene_labels)) 
        gene_labels <- rep("", nrow(y))
    if (!is.null(gmatf$chromosome)) 
        chrom_labels <- gmatf$chromosome
    if (is.null(product_labels)) 
        product_labels <- rep("", nrow(y))
    feature_labels[is.na(feature_labels)] <- ""
    gene_labels[is.na(feature_labels)] <- ""
    product_labels[is.na(feature_labels)] <- ""
    means <- rowMeans(y, na.rm = TRUE)
    means[is.nan(means)] = 0
    y_centered <- y - means
    means2 <- rowMeans(z, na.rm = TRUE)
    means2[is.nan(means2)] = 0
    z_centered <- z - means2
    if (cluster_features == 1) {
        cf <- TRUE
    }
    else {
        cf <- FALSE
    }
    if (feat_span == 1) {
        y_scaled <- y_centered/sqrt(rowMeans(y_centered * y_centered, 
            na.rm = TRUE))
        row_order <- make_o(y_scaled, enable = cf)
    }
    else if (feat_span == 2) {
        z_scaled <- z_centered/sqrt(rowMeans(z_centered * z_centered, 
            na.rm = TRUE))
        row_order <- make_o(z_scaled, enable = cf)
    }
    cluster_samples <- FALSE
    if (clusterby == 1) {
        col_order <- make_o(t(y_centered), enable = cluster_samples)
    }
    else if (clusterby == 2) {
        cluster_samples <- TRUE
        col_order <- make_o(t(y_centered), enable = cluster_samples)
    }
    else if (clusterby == 3) {
        cluster_samples <- TRUE
        col_order <- make_o(t(z_centered), enable = cluster_samples)
    }
    pad <- 0.25
    row_ordering_grob <- ordering_grob(row_order, transpose = TRUE, 
        mirror = TRUE)
    col_ordering_grob <- ordering_grob(col_order)
    heatmap <- heatmap_grob(y_centered[row_order$order, col_order$order, 
        drop = F], signed = TRUE, legend_title = paste0("difference from\nrow mean"))
    heatmap2 <- heatmap_grob(z_centered[row_order$order, col_order$order, 
        drop = F], signed = TRUE, legend_title = paste0("difference from\nrow mean"))
    mean_range <- range(means)
    if (is.na(mean_range[1])) 
        mean_range[1] <- 0.1
    if (is.na(mean_range[2])) 
        mean_range[2] <- 0.1
    if (mean_range[2] == mean_range[1]) 
        mean_range[2] <- mean_range[2] + 1
    mean_range2 <- range(means2)
    if (is.na(mean_range2[1])) 
        mean_range2[1] <- 0.1
    if (is.na(mean_range2[2])) 
        mean_range2[2] <- 0.1
    if (mean_range2[2] == mean_range2[1]) 
        mean_range2[2] <- mean_range2[2] + 1
    mean_graph <- rectGrob(x = rep(mean_range[1], nrow(y)), y = seq_len(nrow(y)) - 
        1, width = means[row_order$order] - mean_range[1], height = rep(1, 
        nrow(y)), just = c(0, 0), default.units = "native", vp = viewport(xscale = mean_range, 
        yscale = c(0, nrow(y))))
    mean_axis <- xaxisGrob(at = axisTicks(mean_range, log = FALSE, 
        nint = 3), label = TRUE, vp = viewport(width = 1, height = 0, 
        y = 1, xscale = mean_range), gp = gpar(cex = 0.75))
    mean_graph2 <- rectGrob(x = rep(mean_range2[1], nrow(z)), 
        y = seq_len(nrow(z)) - 1, width = means2[row_order$order] - 
            mean_range2[1], height = rep(1, nrow(z)), just = c(0, 
            0), default.units = "native", vp = viewport(xscale = mean_range2, 
            yscale = c(0, nrow(z))))
    mean_axis2 <- xaxisGrob(at = axisTicks(mean_range2, log = FALSE, 
        nint = 3), label = TRUE, vp = viewport(width = 1, height = 0, 
        y = 1, xscale = mean_range2), gp = gpar(cex = 0.75))
    mean_label <- textGrob("row\nmean")
    feature_label_grob <- shrinktext_grob(feature_labels[row_order$order], 
        x = rep(0, nrow(y)), y = seq_len(nrow(y)) - 0.5, just = c(0, 
            0.5), vp = viewport(xscale = c(0, 1), yscale = c(0, 
            nrow(y))))
    gene_label_grob <- shrinktext_grob(gene_labels[row_order$order], 
        x = rep(0, nrow(y)), y = seq_len(nrow(y)) - 0.5, just = c(0, 
            0.5), vp = viewport(xscale = c(0, 1), yscale = c(0, 
            nrow(y))))
    product_label_grob <- shrinktext_grob(product_labels[row_order$order], 
        x = rep(0, nrow(y)), y = seq_len(nrow(y)) - 0.5, just = c(0, 
            0.5), vp = viewport(xscale = c(0, 1), yscale = c(0, 
            nrow(y))))
    chrom_label_grob <- shrinktext_grob(chrom_labels[row_order$order], 
        x = rep(0, nrow(y)), y = seq_len(nrow(y)) - 0.5, just = c(0, 
            0.5), vp = viewport(xscale = c(0, 1), yscale = c(0, 
            nrow(y))))
    sample_label_grob <- vertical_shrinktext_grob(sample_labels[col_order$order], 
        x = seq_len(ncol(y)) - 0.5, y = rep(1, ncol(y)), just = c(1, 
            0.5), vp = viewport(xscale = c(0, ncol(y)), yscale = c(0, 
            1)))
    sample_label_grob2 <- vertical_shrinktext_grob(sample_labels2[col_order$order], 
        x = seq_len(ncol(z)) - 0.5, y = rep(1, ncol(z)), just = c(1, 
            0.5), vp = viewport(xscale = c(0, ncol(z)), yscale = c(0, 
            1)))
    frame <- frameGrob(layout = grid.layout(nrow = 4, ncol = 9))
    frame <- packGrob(frame, varistran_grob(col_ordering_grob, 
        height = "inherit", pad = pad), row = 1, col = 2)
    frame <- packGrob(frame, varistran_grob(col_ordering_grob, 
        height = "inherit", pad = pad), row = 1, col = 4)
    lab1 <- textGrob("Mean tail length")
    lab2 <- textGrob("Expression")
    lab3 <- textGrob("Gene Name")
    lab4 <- textGrob("Row Name")
    lab5 <- textGrob("Product Name")
    lab6 <- textGrob("Chr")
    if (cluster_samples == FALSE) {
        frame <- packGrob(frame, varistran_grob(lab1, height = "inherit", 
            pad = pad), row = 1, col = 2)
        frame <- packGrob(frame, varistran_grob(lab2, height = "inherit", 
            pad = pad), row = 1, col = 4)
    }
    frame <- packGrob(frame, varistran_grob(lab3, height = "inherit", 
        pad = pad), row = 1, col = 6)
    frame <- packGrob(frame, varistran_grob(lab4, height = "inherit", 
        pad = pad), row = 1, col = 8)
    frame <- packGrob(frame, varistran_grob(lab5, height = "inherit", 
        pad = pad), row = 1, col = 9)
    frame <- packGrob(frame, varistran_grob(mean_label, height = "inherit", 
        pad = pad), row = 1, col = 3)
    frame <- packGrob(frame, varistran_grob(row_ordering_grob, 
        width = "inherit", pad = pad), row = 2, col = 1)
    frame <- packGrob(frame, varistran_grob(chrom_label_grob, 
        width = "inherit", pad = pad), row = 2, col = 7)
    frame <- packGrob(frame, varistran_grob(lab6, height = "inherit", 
        pad = pad), row = 1, col = 7)
    frame <- packGrob(frame, varistran_grob(heatmap$heatmap, 
        pad = pad), row = 2, col = 2)
    frame <- packGrob(frame, varistran_grob(heatmap2$heatmap, 
        pad = pad), row = 2, col = 4)
    frame <- packGrob(frame, varistran_grob(mean_graph2, width = unit(3, 
        "lines"), pad = pad), row = 2, col = 5)
    frame <- packGrob(frame, varistran_grob(sample_label_grob2, 
        height = "inherit", pad = pad), row = 3, col = 4)
    frame <- packGrob(frame, varistran_grob(mean_label, height = "inherit", 
        pad = pad), row = 1, col = 5)
    frame <- packGrob(frame, varistran_grob(mean_axis2, height = unit(3, 
        "lines"), pad = pad), row = 3, col = 5)
    frame <- packGrob(frame, varistran_grob(mean_graph, width = unit(3, 
        "lines"), pad = pad), row = 2, col = 3)
    frame <- packGrob(frame, varistran_grob(feature_label_grob, 
        width = "inherit", pad = pad), row = 2, col = 8)
    frame <- packGrob(frame, varistran_grob(gene_label_grob, 
        width = "inherit", pad = pad), row = 2, col = 6)
    frame <- packGrob(frame, varistran_grob(product_label_grob, 
        width = "inherit", pad = pad), row = 2, col = 9)
    frame <- packGrob(frame, varistran_grob(sample_label_grob, 
        height = "inherit", pad = pad), row = 3, col = 2)
    frame <- packGrob(frame, varistran_grob(mean_axis, height = unit(3, 
        "lines"), pad = pad), row = 3, col = 3)
    frame <- packGrob(frame, varistran_grob(heatmap$legend, height = "inherit", 
        pad = pad), row = 4, col = 3)
    frame <- packGrob(frame, varistran_grob(heatmap2$legend, 
        height = "inherit", pad = pad), row = 4, col = 5)
    outer <- frameGrob()
    outer <- packGrob(outer, varistran_grob(frame), row = 1, 
        col = 1)
    result <- varistran_grob(outer, pad = pad)
    result$info <- list(row_order = row_order, col_order = col_order)
    result
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
