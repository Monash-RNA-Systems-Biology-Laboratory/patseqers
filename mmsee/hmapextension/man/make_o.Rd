\name{make_o}
\alias{make_o}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make_o(mat, enable = TRUE, fast = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
%%     ~~Describe \code{mat} here~~
}
  \item{enable}{
%%     ~~Describe \code{enable} here~~
}
  \item{fast}{
%%     ~~Describe \code{fast} here~~
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
function (mat, enable = TRUE, fast = FALSE) 
{
    mat[is.na(mat)] = 0.001
    if (nrow(mat) < 3 || !enable) {
        list(dendrogram = NULL, order = seq_len(nrow(mat)), paths = rep("", 
            nrow(mat)))
    }
    else if (fast) {
        perm <- seriation::seriate(mat, method = "PCA")[[1]]
        list(dendrogram = NULL, order = seriation::get_order(perm), 
            paths = rep("", nrow(mat)))
    }
    else {
        dist_mat <- dist(mat)
        control <- list(hclust = hclust(dist_mat))
        dend_mat <- as.dendrogram(seriation::seriate(dist_mat, 
            method = "OLO", control = control)[[1]])
        list(dendrogram = dend_mat, order = order.dendrogram(dend_mat), 
            paths = dendrogram_paths(dend_mat))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
