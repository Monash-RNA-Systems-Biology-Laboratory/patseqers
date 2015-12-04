#' @export
sh_hmap_detailed <-
    function(rw, sample_labels=NULL, sample_labels2=NULL, feature_labels=NULL, prefix="") {
        p <- function(name) paste0(prefix,name)
        sample_labels <- ensure_reactable(sample_labels)
        sample_labels2 <- ensure_reactable(sample_labels2)
        feature_labels <- ensure_reactable(feature_labels)
        plot <- shiny_plot(
            callback = function(env) {
                print(env[[p("grob")]]())
            },
            width=1250,
            height=900,
            dlname="heatmap",
            prefix=p("plot_")
        )
        
        ui <- shiny::tags$div(
            shiny::titlePanel("Heatmap"),
            shiny::p("Features are selected based on span of:"),
            shiny::radioButtons(p("featspan"), 
                                label="Expression or Tail length", 
                                choices=list("Tail Length"=1, "Expression"=2), 
                                selected=1,
                                inline=TRUE),
            shiny::uiOutput(p("chrs")),
            shiny::numericInput(p("n"), "Number of features to show", 50, min=10,max=2000,step=10),
            shiny::radioButtons(p("selFeat"), 
                                label="Cluster by:", 
                                choices=list("Tail count"=1, "Expression"=2), 
                                selected=1,
                                inline=TRUE),
            shiny::numericInput(p("nmin"), "Trim low Tail_count to NA", 10, min=0,max=1000,step=1),
            shiny::radioButtons(p("seqGroup"), 
                                label="Cluster/Group options", 
                                choices=list("Cluster by row"=1, "Group by location"=2), 
                                selected=1,
                                inline=TRUE),
            shiny::radioButtons("Clusterby", 
                                label="Cluster by samples: ", 
                                choices=list("None" = 1, "Tail length" = 2, "Expression" = 3), 
                                selected = 1,
                                inline=TRUE),
            plot$component_ui,
            parenthetically("This plot is produced by a modified varistran::plot_heatmap.")
        )
        #Processes the input list into a single dataframe with annotation and count data
        wproc <- function(env){
            proc <- reactive({
                rw2 <- list()
                
                rw$Tail[rw$Tail_count < env$input[[p("nmin")]]] = NA
                rw2$Tail <- rw$Tail
                rw2$Count <- rw$Count
                #Append names to grep out
                colnames(rw2$Tail) <- paste(colnames(rw2$Tail), "Tail")
                colnames(rw2$Count) <- paste(colnames(rw2$Count), "Count")
                rw3 <- list()
                
                rw3$Tail <- rw2$Tail
                rw3$Count <- varistran::vst(rw2$Count)
                
                rwm <- list()
                rwm$data <- merge(rw3$Tail, rw3$Count, by=0)
                rownames(rwm$data) <- rwm$data$Row.names
                rw3$annotate<- rw$Annotation[rownames(rw$Annotation) %in% rownames(rw3$Tail),]
                rw3$Tail_count <- rw$Tail_count[rownames(rw$Tail_count) %in% rownames(rw3$Tail),]
                rw3$annotate <- rw3$annotate[order(rw3$annotate$chromosome, rw3$annotate$start),]
                rw3$Tail <- rw3$Tail[order(rw3$annotate$chromosome, rw3$annotate$start),]
                rw3$Count <- rw3$Count[order(rw3$annotate$chromosome, rw3$annotate$start),]
                rw3$Tail_count <- rw3$Tail_count[order(rw3$annotate$chromosome, rw3$annotate$start),]
                #Sort into order
                cvec <- rw3$annotate$chromosome
                cvec <- cvec %in% env$input[[p("choosechr")]]
                
                rw3$annotate <- rw3$annotate[cvec,]
                rw3$Tail <- rw3$Tail[cvec,]
                rw3$Tail_count <- rw3$Tail_count[cvec,]
                rw3$Count <- rw3$Count[cvec,]
                rownames(rw3$Count) <- substr(rownames(rw3$Count), 1, 17)
                rownames(rw3$Tail) <- substr(rownames(rw3$Tail), 1, 17)
                rownames(rw3$Tail_count) <- substr(rownames(rw3$Tail_count), 1, 17)
                rownames(rw3$annotate) <- substr(rownames(rw3$annotate), 1, 17)
                rw3$annotate$product <- substr(rw3$annotate$product , 1, 76)
                rw3$annotate$gene <- substr(rw3$annotate$gene, 1, 11)
                #Truncate names for neatness
                return(rw3)
            })
            return(proc())
        }
        server <- function(env) {
            
            env$output[[p("chrs")]] <- shiny::renderUI({
                tmpvec <- levels(rw$Annotation$chromosome)
                selectizeInput("choosechr", "Choose chromosomes to display",multiple=T,tmpvec, selected=tmpvec)
            })
            env[[p("grob")]] <- reactive({
                
                a1 <- wproc(env)
                if(env$input[[p("selFeat")]] == 1){
                    y <- a1$Tail
                } else if(env$input[[p("selFeat")]] == 2){
                    y <- a1$Count
                }
                
                y <- ensure_reactable(y)
                
                n <- env$input[[p("n")]]
                if (n > 2000) stop("Drawing large heatmaps uses excessive system resources. Sorry.")
                
                y_val <- as.matrix(y(env))
                y_span <- apply(y_val,1,max) - apply(y_val,1,min)
                selection <- rep(FALSE,nrow(y_val))
                selection[ order(-y_span)[ seq_len(n) ] ] <- TRUE
                
                if (sum(selection) < 1) stop("No features to show.")
                
                pl_hmap_detailed(
                    matf1=a1$Tail[selection,,drop=FALSE],
                    matf2=a1$Count[selection,,drop=FALSE],
                    gmatf=a1$annotate[selection,,drop=FALSE],
                    sample_labels=sample_labels(env),
                    sample_labels2=sample_labels(env),
                    feature_labels=feature_labels(env)[selection],
                    clusterby=env$input[[p("Clusterby")]],
                    cluster_features=env$input[[p("seqGroup")]],
                    feat_span=env$input[[p("featspan")]]
                )
            })
            
            plot$component_server(env)
        }
        
        composable_shiny_app(ui, server)
    }
