#' @title Produces detailed heatmap
#' @details 
#' Workhorse function for this package.
#'      
#' @param rw List of dataframes 
#' Takes a read.grouped.table() as input or a list of four dataframes (more data frames are ok but it only uses these):
#' Counts - Genewise counts of expression
#' Tail - Mean tail length
#' Tail_counts - Number of poly-A tails counted
#' Annotation - Information regarding the annotation information
#'      Gene name, chromosome, gene product, biotype etc...
#'  
#' @param sample_labels Sample labels
#' @param sample_labels2 Sample labels (second plot)
#' @param feature_labels Feature labels
#' @param prefix Prefix for plot
#' @param species Species of the data. Currently supports Human (Hs), Saccharomyces cerevisiae (Sc), Caenorhabditis elegans (Ce), Mus musculus (Mm)
#' 
#' @return 
#' Returns a composable shiny app object
#' 
#' @import shiny
#' @import varistran
#' @import shinyURL
#' @export

sh_hmap_detailed <- function(rw, sample_labels=NULL, sample_labels2=NULL, feature_labels=NULL, prefix="", gotable=NULL, species=NULL) {
    if(is.null(species) || sum(species == c("Hs", "Sc", "Ce", "Mm"))!=1)
        stop("Species argument missing or incorrect")
    
    p <- function(name) paste0(prefix,name)
    sample_labels <- ensure_reactable(sample_labels)
    sample_labels2 <- ensure_reactable(sample_labels2)
    feature_labels <- ensure_reactable(feature_labels)
    gotable <- ensure_reactable(gotable)
    
    plot <- shiny_p(
        callback = function(env) {
            print(env[[p("grob")]]())
        },
        rorder = function(env) {
            env[[p("grob")]]()$info$row_order$order
        },
        width=1250,
        height=900,
        dlname="heatmap",
        prefix=p("plot_"),
        selin = function(env){
            env$seldat()$ann
        },
        spp=species
    )
    # Shiny's UI layout 
    ui <- shiny::tags$div(
        shiny::titlePanel("Heatmap"),
        shiny::tabsetPanel(
            shiny::tabPanel("Options",
                            shiny::fluidRow(
                                shiny::column(3,
                                              shiny::p("Features are selected based on span of:"),
                                              shiny::radioButtons(p("selFeat"), 
                                                                  label="Select features by:", 
                                                                  choices=list("Tail count"=1, "Expression"=2), 
                                                                  selected=1,
                                                                  inline=TRUE),
                                              shiny::uiOutput(p("chrs"))
                                ),
                                shiny::column(3,
                                              shiny::numericInput(p("n"), "Number of features to show", 50, min=10,max=2000,step=10),
                                              shiny::numericInput(p("nmin"), "Trim Tail Counts below value to NA", 5, min=0,max=1000,step=1),
                                              shiny::numericInput(p("expmin"), "Exclude rows with low expression counts", 0, min=0,max=1500,step=1),
                                              shiny::radioButtons(p("roword"), 
                                                                  label="Features ordered by: ", 
                                                                  choices=list("Tail Length"=1, "Expression"=2, "Group by location"=3), 
                                                                  selected=1,
                                                                  inline=TRUE),
                                              shiny::radioButtons(p("clusterby"), 
                                                                  label="Order samples by: ", 
                                                                  choices=list("Tail length" = 2, "Expression" = 3, "Manually order from select columns" = 1), 
                                                                  selected = 1,
                                                                  inline=TRUE)
                                ),
                                shiny::column(3,
                                              shiny::uiOutput(p("selCol")),
                                              shinyURL::shinyURL.ui()
                                )
                                
                            )),
            shiny::tabPanel("Plot",plot$component_ui),
            shiny::tabPanel("GO analysis of selection",
                            shiny::tableOutput("tabout"))
        ),
        parenthetically("This plot is produced by a modified varistran::plot_heatmap.")
        
    )
    
    # Shiny's server
    server <- function(env) {
        shinyURL::shinyURL.server(env$session)
        # Processes the input from rw into the correct length and returns a list of 4 data frames
        wproc <- reactive({
            
            rw2 <- list()
            
            colvec <- -which(names(rw$Tail) %in% env$input[[p("choosecol")]])
            
            rw$Tail[rw$Tail_count < env$input[[p("nmin")]]] = NA
            
            rw2$Tail <- rw$Tail[,-colvec]
            rw2$Count <- rw$Count[,-colvec]
            
            hldvec2 <- apply(rw2$Tail,1,function(row) {
                any(sum(!is.na(row))>= 2)
            })
            
            rw2$Tail <- rw2$Tail[hldvec2,]
            rw2$Count <- rw2$Count[hldvec2,]
            
            if(env$input[[p("clusterby")]] == 4){
                rw$Tail <- rw$Tail[env$input[[p("choosecol")]]]
                rw2$Count <- rw$Count[env$input[[p("choosecol")]]]
            }
            
            #Append names for neatness in graphic
            colnames(rw2$Tail) <- paste(colnames(rw2$Tail), "Tail")
            colnames(rw2$Count) <- paste(colnames(rw2$Count), "Count")
            rw3 <- list()
            
            rw3$Tail <- rw2$Tail
            rw3$Count <- varistran::vst(rw2$Count)
            
            rw3$annotate<- rw$Annotation[rownames(rw$Annotation) %in% rownames(rw3$Tail),]
            
            rw3$Tail_count <- rw$Tail_count[rownames(rw$Tail_count) %in% rownames(rw3$Tail),]
            
            orderedvec <- order(rw3$annotate$chromosome, rw3$annotate$start)
            
            rw3$annotate <- rw3$annotate[orderedvec,]
            rw3$Tail <- rw3$Tail[orderedvec,]
            rw3$Count <- rw3$Count[orderedvec,]
            rw3$Tail_count <- rw3$Tail_count[orderedvec,]
            rw3$Tail_count <- rw3$Tail_count[,-colvec]
            
            #Sort into order
            cvec <- rw3$annotate$chromosome
            cvec <- cvec %in% env$input[[p("choosechr")]]
            
            rw3$annotate <- rw3$annotate[cvec,]
            rw3$Tail <- rw3$Tail[cvec,]
            rw3$Tail_count <- rw3$Tail_count[cvec,]
            rw3$Count <- rw3$Count[cvec,]
            
            #Truncate names for neatness
            rownames(rw3$Count) <- substr(rownames(rw3$Count), 1, 17)
            rownames(rw3$Tail) <- substr(rownames(rw3$Tail), 1, 17)
            rownames(rw3$Tail_count) <- substr(rownames(rw3$Tail_count), 1, 17)
            rownames(rw3$annotate) <- substr(rownames(rw3$annotate), 1, 17)
            rw3$annotate$product <- substr(rw3$annotate$product , 1, 76)
            rw3$annotate$gene <- substr(rw3$annotate$gene, 1, 11)
            
            #Trim out values if we want to
            hldvec2 <- apply(rw3$Count,1,function(row) {
                any(row >= log2(env$input[[p("expmin")]]))
            })
            
            rw3$annotate <-data.frame(rw3$annotate[hldvec2,])
            rw3$Tail_count <- data.frame(rw3$Tail_count[hldvec2,])
            rw3$Tail <- data.frame(rw3$Tail[hldvec2,])
            rw3$Count <- data.frame(rw3$Count[hldvec2,])
            
            return(rw3)
        })
        
        # Processes the selection of rows by calculating maximum span in either tail length of expression
        selproc <- reactive({
            a1 <- wproc()
            if(env$input[[p("selFeat")]] == 1){
                y <- a1$Tail
            } else if(env$input[[p("selFeat")]] == 2){
                y <- a1$Count
            }
            y <- ensure_reactable(y)
            
            n <- env$input[[p("n")]]
            if (n > 2000) stop("Drawing large heatmaps uses excessive system resources. Sorry.")
            
            y_val <- as.matrix(y(env))
            y_span <- apply(y_val,1,max,na.rm=T) - apply(y_val,1,min,na.rm=T)
            selection <- rep(FALSE,nrow(y_val))
            selection[ order(-y_span)[ seq_len(n) ] ] <- TRUE
            
            if (sum(selection) < 1) stop("No features to show.")
            
            rtVal <- list()
            rtVal$sel <- selection
            rtVal$ann <- a1$annotate[selection,,drop=FALSE]
            rtVal$a1 <- a1
            return(rtVal)
        })        
        
        goDBret <- reactive({
            return(env$output[[p("GOsrch")]])
        })
        # Enables env to hold the data from sh_hmap_detailed
        env$seldat <- selproc
        env$goDB <- goDBret
        
        env$output[[p("debugOut")]] <- shiny::renderText({
            names(env)
        })
        # RenderUI output for shiny, dynamically generate two selctize elements ---
        env$output[[p("chrs")]] <- shiny::renderUI({
            tmpvec <- levels(rw$Annotation$chromosome)
            selectizeInput("choosechr", "Choose chromosomes to display",multiple=T,tmpvec, selected=tmpvec)
        })
        env$output[[p("selCol")]] <- shiny::renderUI({
            colvec <- names(rw$Tail)
            selectizeInput("choosecol", "Choose samples to display",multiple=T,colvec, selected=colvec)
        })
        #---
        
        env[[p("grob")]] <- reactive({
            
            selrt <- selproc()
            splitDF <- selrt$a1
            selection <- selrt$sel
            
            pl_hmap_detailed(
                matf1=splitDF$Tail[selection,,drop=FALSE],
                matf2=splitDF$Count[selection,,drop=FALSE],
                gmatf=splitDF$annotate[selection,,drop=FALSE],
                sample_labels=sample_labels(env),
                sample_labels2=sample_labels(env),
                feature_labels=feature_labels(env)[selection],
                clusterby=env$input[[p("clusterby")]],
                row_ord=env$input[[p("roword")]]   
            )
        })
        plot$component_server(env)
        
        env$output[[p("tabout")]] <- shiny::renderTable({
            env$gotab()
        })
    }
    composable_shiny_app(ui, server)
}
