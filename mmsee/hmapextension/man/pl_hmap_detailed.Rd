% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pl_hmap_detailed.R
\name{pl_hmap_detailed}
\alias{pl_hmap_detailed}
\title{Printable heatmap grob}
\usage{
pl_hmap_detailed(matf1, matf2, gmatf, clusterby = 1, sample_labels = NULL,
  sample_labels2 = NULL, feature_labels = NULL, gene_labels = NULL,
  product_labels = NULL, row_ord = 1)
}
\arguments{
\item{matf1}{Data frame of Tail length}

\item{matf2}{Data frame of Counts (genewise expression) (Should already be normalised)}

\item{gmatf}{Data frame of annotation data}

\item{clusterby}{Cluster columns by tail length or expression (defaults to None)}

\item{sample_labels}{Sample labels}

\item{sample_labels2}{Sample labels (second plot)}

\item{feature_labels}{Feature labels}

\item{gene_labels}{Gene labels}

\item{product_labels}{Product labels}

\item{row_ord}{Order rows by tail length or expression}
}
\details{
Produces a print-able heatmap grob
}

