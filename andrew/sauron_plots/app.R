library("shiny")
library("nesoni")
library("ggplot2")
library("reshape2")

ui <- fluidPage (
  titlePanel("Sauron Plotter"),
  
  em(helpText("created by Andrew Pattison and Paul Harrison for the Beilharz Lab", align = "right")),
  
  helpText("This app plots gene expressions vs tail length and peak-pair expresion shift"),
  
  br(),
  
  sidebarPanel(    
    uiOutput("select_file_path"),
    
    selectInput("select_plot_meth", label = h5("What do you want to plot?"), 
                choices = list("Counts vs Tail Length" = 1, "Counts vs Peak Shift" = 2, 
                               "Peak Shift vs Tail Length" = 3), selected = 1),
    
    checkboxInput("tmm_norm", label = "TMM Normalise Counts", value = F),
    
    uiOutput("choose_samples")
    
  ),
  
  mainPanel(
    tabsetPanel(
      tabPanel("Plot",
               h2("Plot Output",  height = "600"),   
               plotOutput('plot', brush = "plot_brush")
      ),
      tabPanel("Plot Frame",
               h2("The data frame that the plots are made from",  height = "600"),   
               dataTableOutput("data")
      ),
      tabPanel("Raw Count Data",
               h2("The raw counts data frame",  height = "600"),   
               dataTableOutput("raw_count_data")
      
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
  
  output$raw_count_data <- renderDataTable({
    count_info_table()    
  })
  
  c_v_t_info_table<- reactive({
    get_plot_cols_c_v_t (count_info_table(), input$tmm_norm)
  })
  
  samples_from_counts<- reactive({
    get_sample_names (c_v_t_info_table())
  })
  
  output$choose_samples <- renderUI({
    checkboxGroupInput("select_samples", 
    label = h5("Select which samples you would like to plot"),
                       choices = samples_from_counts(), 
                       selected = samples_from_counts()[1]) 
  })
  
  output$data <- renderDataTable({
    if (input$select_plot_meth == 1){
      return(c_v_t_info_table())    
    }
    else if (input$select_plot_meth == 2){
      return(count_info_table()) 
    }
    else{
      return(count_info_table())
    }
  })
  
  plot_calcs <- reactive({
    if (input$select_plot_meth == 1){
      plot <- make_plot_c_v_t(c_v_t_info_table(),input$select_samples)
      return(plot)    
    }
    else if (input$select_plot_meth == 2){
      return(count_info_table()) 
    }
    else{
      return(count_info_table())
    }    
  }) 
  output$plot <- renderPlot({  
    plot_calcs()
  })
  
}

make_info_frame <- function (tt_location){
  counts_csv <- read.grouped.table(paste0(tt_location, "/counts.csv"))  
  return (data.frame(counts_csv))
}

get_plot_cols_c_v_t <- function (count_table, tmm_norm){
  if (tmm_norm == F){
    counts <- count_table [,grep("Count.*",colnames(count_table))]
    mean_tails <- count_table [,grep("Tail.*",colnames(count_table))]
    final <- data.frame(count_table$Annotation.gene, counts, mean_tails)
    return(final)
  }
  
}

get_sample_names <- function(df){
  
  unrefined_names <- colnames(df[,grep("Count.*", colnames(df))])
  refined_names <- list()
  for (name in unrefined_names){
    refined_name <- strsplit(name, "[.]") [[1]][[2]]
    refined_names <- c(refined_names, refined_name)
  }
  return(refined_names)
}

make_plot_c_v_t <- function (df, samples){
  count_frame <- data.frame()
  for (sample in samples){
    sample_col_count <- df[,grep(paste0("Count.", sample), colnames(df))] 
    sample_col_tail <- df[,grep(paste0("Tail.", sample), colnames(df))] 
    to_bind <- data.frame(sample_col_count ,sample_col_tail, sample)
    str(to_bind)
    count_frame <- rbind(count_frame,to_bind)

  }  
  return(
    ggplot(data=count_frame, aes (x=log2(sample_col_count), y=sample_col_tail))+
           geom_point(aes(colour= sample, group= sample))
    )
}

shinyApp (ui = ui, server = server)
