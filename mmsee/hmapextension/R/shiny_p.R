#' @title Integrates heatmap into shiny
#' @import gridBase

shiny_p <- function(callback, width=500, height=500, dlname="plot", prefix="", selin, rorder) {
    selin <- ensure_reactable(selin)
    rorder <- ensure_reactable(rorder)
    
    p <- function(name) paste0(prefix,name)
    ui <- shiny::tags$div(
        shiny::fluidRow(
            shiny::column(3, shiny::numericInput(p("width"), "Plot width", width, min=100, max=10000, step=50)),
            shiny::column(3, shiny::numericInput(p("height"), "Plot height", height, min=100, max=10000, step=50)),
            shiny::column(4, shiny::tags$label("Download"), shiny::tags$br(),
                          shiny::downloadButton(p("pdf"), "PDF"),
                          shiny::downloadButton(p("eps"), "EPS"))
        ),
        shiny::uiOutput(p("plotui"), width="auto", height="auto"),
        DT::dataTableOutput(p('datab'))
    )
    
    server <- function(env) {
        output <- env$output
        
        i <- function(name) env$input
        
        output[[p("plotui")]] <- shiny::renderUI({
            
            shiny::plotOutput(p("plot"),
                              brush = brushOpts(
                                  id ="plot_brush",  
                                  direction = 'y',
                                  resetOnNew = T,
                                  clip=T
                              ),
                              width=env$input[[p("width")]],
                              height=env$input[[p("height")]]
            )
        })
        
        calcdt <- reactive({
            numrows <- nrow(selin(env))
            
            ytop <- round((env$input$plot_brush$ymax + 2)/(54/numrows))
            ybot <- round((env$input$plot_brush$ymin + 2)/(54/numrows)+1)
            if(length(ytop) == 0){
                return(selin(env)[rorder(env),])
            } else {
                sel <- selin(env)[rorder(env),]
                sel <- sel[(ytop):ybot,]
                return(sel)
            }
        })
            
        output[[p("datab")]] <- DT::renderDataTable(calcdt(), server=F,
                                                        options = list(searchHighlight = TRUE)
        )

        
        output[[p("plot")]] <- shiny::renderPlot({ 
            vp <- viewport(layout.pos.row = 1, layout.pos.col = 1)
            plot.new()
            callback(env)
            seekViewport("prodVP")
            par(new=T, plt=gridPLT())
            plot(1,type="n", axes=F, xlab ="",ylab="",xlim=c(0,50),ylim=c(0,50))
            
            popViewport()
            
        }
        )
        output[[p("pdf")]] <- shiny::downloadHandler(
            paste0(dlname,".pdf"),
            function(filename) {
                pdf(filename, width=i("width")/72, height=i("height")/72)
                callback(env)
                dev.off()
            }
        )
        output[[p("eps")]] <- shiny::downloadHandler(
            paste0(dlname,".eps"),
            function(filename) {
                postscript(filename, width=i("width")/72, height=i("height")/72,
                           paper="special", onefile=FALSE, horizontal=FALSE)
                callback(env)
                dev.off()
            }
        )
    }
    
    composable_shiny_app(ui, server)
}