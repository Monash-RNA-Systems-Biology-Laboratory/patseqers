library("shiny")
library("nesoni")
library("ggplot2")
library("reshape2")
library("varistran")

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
    
    checkboxInput("varis", label = "Varistran Tansform Counts", value = T),
    checkboxInput("combine", label = "Combine by Replicates", value = F),
    
    uiOutput("choose_samples"),
    uiOutput("select_group")
    
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
    get_plot_cols_c_v_t (count_info_table(), input$varis, input$combine)
  })
  
  samples_from_counts<- reactive({
    get_sample_names (count_info_table())
  })
  select_group_fun <- reactive({
    count <- 1
    lapply(input$select_samples, 
           function(i) {          
             selectInput(paste0('snumber', i),              
                         h5(paste0('Select a group for ', i)),
                         choices = 1:length(input$select_bam_files))
           }
    )
  })  
  group_list <- reactive({
    res <- lapply(input$select_samples, 
                  function(i) { 
                    input[[paste0('snumber', i)]]
                  })
  })
  output$select_group <- renderUI({
    select_group_fun()
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
      plot <- make_plot_c_v_t(c_v_t_info_table(),input$select_samples, 
                              input$select_plot_meth, input$file_path)
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

get_plot_cols_c_v_t <- function (count_table, varis, combine){
  print(head(count_table))
  if (combine == T){
    
  }
  if (varis == F){
    counts <- count_table [,grep("Count.*",colnames(count_table))]

  }
  else{
    counts <- varistran::vst(count_table [,grep("Count.*",colnames(count_table))])
    
  }
  mean_tails <- count_table [,grep("Tail.*",colnames(count_table))]
  final <- data.frame(count_table$Annotation.gene, counts, mean_tails)
  return(final)
  
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

make_plot_c_v_t <- function (df, samples, name, title){
  if (name ==1){
    t <- paste ("Counts vs Tail Length",title)
    x <- "Count"
    y <- "Tail Length"
  }
  else if (name ==2){
    t <- paste ("Counts vs Peak Shift",title)
    x <- "Count"
    y <- "Peak Shift"
  }
  else{
    t <- paste ("Peak Shift vs Tail Length",title)
    x <- "Peak Shift"
    y <- "Tail Length"    
  }
  count_frame <- data.frame()
  for (sample in samples){
    sample_col_count <- df[,grep(paste0("Count.", sample), colnames(df))] 
    sample_col_tail <- df[,grep(paste0("Tail.", sample), colnames(df))] 
    to_bind <- data.frame(sample_col_count ,sample_col_tail, sample)
    count_frame <- rbind(count_frame,to_bind)

  }  
  return(
    ggplot(data=count_frame, aes (x=sample_col_count, y=sample_col_tail))+
           geom_point(aes(colour= sample, group= sample))+
            labs(title = t, x = x, y = y)
    )
}

shinyApp (ui = ui, server = server)
