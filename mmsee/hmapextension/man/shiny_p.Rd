% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shiny_p.R
\name{shiny_p}
\alias{shiny_p}
\title{Integrates heatmap into shiny}
\usage{
shiny_p(callback, width = 500, height = 500, dlname = "plot",
  prefix = "", selin, rorder, spp, hgc = 0.05, otype = 1,
  goenabl = TRUE)
}
\arguments{
\item{callback}{Heatmap Grob from pl_hmap_detailed}

\item{width}{Default width of the heatmap grob when the shiny app loads}

\item{height}{Default height of the heatmap grob when the shiny app loads}

\item{dlname}{Default filename when downloading the heatmap image}

\item{prefix}{Prefix}

\item{selin}{Selection of genes from sh_hmap_detailed}

\item{rorder}{Row order from sh_hmap_detailed}

\item{spp}{Species of the data (Hs, Sc, Ce, Mm)}

\item{goenabl}{Whether GO Term analysis is enabled. TRUE by default. This might need changing in future}
}
\details{
Shiny wrapper for sh_hmap_detailed
Uses gridBase to produce brush for interactivity with plot
}

