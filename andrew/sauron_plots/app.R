library("shiny")
library("nesoni")

ui <- fluidPage (
  titlePanel("Sauron Plotter"),
  
  em(helpText("created by Andrew Pattison and Paul Harrison for the Beilharz Lab", align = "right")),
  
  helpText("This app plots gene expressions vs tail length and peak-pair expresion shift"),
  
  br(),
  
  sidebarPanel(    
    uiOutput("select_file_path"),
    
    selectInput("select_plot_meth", label = h3("What do you want to plot?"), 
                choices = list("Counts vs Tail Length" = 1, "Counts vs Peak Shift" = 2, 
                               "Peak Shift vs Tail Length" = 3), selected = 1)
  ),
  
  mainPanel(
    tabsetPanel(
      tabPanel("Plot",
               h2("Plot Output",  height = "600"),   
               plotOutput('scp_plot', brush = "plot_brush")
      ),
      tabPanel("Info",
               h2("Info",  height = "600"),   
               dataTableOutput("data")
      )
    )
  )
)

server <- function (input, output){
  
  output$select_file_path <- renderUI({
    selectInput("file_path", label = h4("Select a Tail-Tools Output"), 
                choices = list.dirs(full.names=F, recursive =F), 
                selected =  list.dirs(full.names=F, recursive =F)[1])   
  })
  
  count_info_table <- reactive({
    make_info_frame(input$file_path)      
  })
  
  c_v_t_info_table<- reactive({
    get_plot_cols_c_v_t (count_info_table())
  })
  
  output$data <- renderDataTable({
    if (input$select_plot_meth == "Counts vs Tail Length"){
      return(c_v_t_info_table())    
    }
    else if (input$select_plot_meth == "Counts vs Peak Shift"){
      return(count_info_table()) 
    }
    else{
      return(count_info_table())
    }
  })
}

make_info_frame <- function (tt_location){
  counts_csv <- read.grouped.table(paste0(tt_location, "/counts.csv"))  
  return (data.frame(counts_csv))
}

get_plot_cols_c_v_t <- function (count_table){
  
}

shinyApp (ui = ui, server = server)


# make_required_frame <- function (counts_csv, ratios_csv){
#   counts_csv <- data.frame(counts_csv)
#   ratios_csv <- data.frame(ratios_csv)
#   combined_frame <- 
#   return(combined_frame)
# }
# ratios_csv <- read.grouped.table(paste0(tt_location, "/individual-pairs.csv"))
# required_data_frame <- make_required_frame(counts_csv, ratios_csv)
