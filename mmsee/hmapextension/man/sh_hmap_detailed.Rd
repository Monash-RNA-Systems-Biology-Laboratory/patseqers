\name{sh_hmap_detailed}
\alias{sh_hmap_detailed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{pl_hmap_detailed
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sh_hmap_detailed(rw, sample_labels = NULL, sample_labels2 = NULL, feature_labels = NULL, prefix = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rw}{
%%     ~~Describe \code{rw} here~~
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
  \item{prefix}{
%%     ~~Describe \code{prefix} here~~
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
function (rw, sample_labels = NULL, sample_labels2 = NULL, feature_labels = NULL, 
    prefix = "") 
{
    p <- function(name) paste0(prefix, name)
    sample_labels <- ensure_reactable(sample_labels)
    sample_labels2 <- ensure_reactable(sample_labels2)
    feature_labels <- ensure_reactable(feature_labels)
    plot <- shiny_plot(callback = function(env) {
        print(env[[p("grob")]]())
    }, width = 1250, height = 900, dlname = "heatmap", prefix = p("plot_"))
    ui <- shiny::tags$div(shiny::titlePanel("Heatmap"), shiny::p("Features are selected based on span of:"), 
        shiny::radioButtons(p("featspan"), label = "Expression or Tail length", 
            choices = list(Expression = 2, `Tail Length` = 1), 
            selected = 1, inline = TRUE), shiny::uiOutput(p("chrs")), 
        shiny::numericInput(p("n"), "Number of features to show", 
            50, min = 10, max = 2000, step = 10), shiny::numericInput(p("nmin"), 
            "Trim low Tail_count to NA", 10, min = 0, max = 1000, 
            step = 1), shiny::radioButtons(p("seqGroup"), label = "Cluster/Group options", 
            choices = list(`Cluster by row` = 1, `Group by location` = 2), 
            selected = 1, inline = TRUE), shiny::radioButtons("Clusterby", 
            label = "Cluster by samples: ", choices = list(None = 1, 
                `Tail length` = 2, Expression = 3), selected = 1, 
            inline = TRUE), plot$component_ui, parenthetically("This plot is produced by a modified varistran::plot_heatmap."))
    wproc <- function(env) {
        proc <- reactive({
            rw2 <- list()
            rw$Tail[rw$Tail_count < env$input[[p("nmin")]]] = NA
            rw2$Tail <- rw$Tail
            rw2$Count <- rw$Count
            colnames(rw2$Tail) <- paste(colnames(rw2$Tail), "Tail")
            colnames(rw2$Count) <- paste(colnames(rw2$Count), 
                "Count")
            rw3 <- list()
            rw3$Tail <- rw2$Tail
            rw3$Count <- varistran::vst(rw2$Count)
            rwm <- list()
            rwm$data <- merge(rw3$Tail, rw3$Count, by = 0)
            rownames(rwm$data) <- rwm$data$Row.names
            rwm$annotate <- rw$Annotation[rownames(rw$Annotation) \%in\% 
                rwm$data$Row.names, ]
            ra <- rw$Tail_count[rownames(rw$Tail_count) \%in\% 
                rwm$data$Row.names, ]
            colnames(ra) <- paste(colnames(ra), "tc")
            rwm$data$Row.names <- NULL
            hld <- merge(rwm$data, rwm$annotate, by = 0)
            hld <- hld[order(hld$chromosome, hld$start), ]
            cvec <- hld[, grepl("chromosome", colnames(hld))]
            cvec <- cvec \%in\% env$input[[p("choosechr")]]
            hld <- hld[cvec, ]
            rownames(hld) <- hld$Row.names
            rownames(hld) <- substr(rownames(hld), 1, 17)
            hld$product <- substr(hld$product, 1, 76)
            hld$gene <- substr(hld$gene, 1, 11)
            hld$Row.names <- NULL
            return(hld)
        })
        return(proc())
    }
    server <- function(env) {
        env$output[[p("chrs")]] <- shiny::renderUI({
            tmpvec <- levels(rw$Annotation$chromosome)
            selectizeInput("choosechr", "Choose chromosomes to display", 
                multiple = T, tmpvec, selected = tmpvec)
        })
        env[[p("grob")]] <- reactive({
            a1 <- wproc(env)
            tvec <- grepl("Tail|Count", colnames(a1))
            y <- ensure_reactable(a1[, tvec])
            n <- env$input[[p("n")]]
            if (n > 2000) 
                stop("Drawing large heatmaps uses excessive system resources. Sorry.")
            y_val <- as.matrix(y(env))
            y_span <- apply(y_val, 1, max) - apply(y_val, 1, 
                min)
            selection <- rep(FALSE, nrow(y_val))
            selection[order(-y_span)[seq_len(n)]] <- TRUE
            if (sum(selection) < 1) 
                stop("No features to show.")
            pl_hmap_detailed(matf = y_val[selection, , drop = FALSE], 
                gmatf = a1[!tvec][selection, , drop = FALSE], 
                sample_labels = sample_labels(env), sample_labels2 = sample_labels(env), 
                feature_labels = feature_labels(env)[selection], 
                clusterby = env$input[[p("Clusterby")]], cluster_features = env$input[[p("seqGroup")]], 
                feat_span = env$input[[p("featspan")]])
        })
        plot$component_server(env)
    }
    composable_shiny_app(ui, server)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
