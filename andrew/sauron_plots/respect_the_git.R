library("shiny")
library("nesoni")
library("ggplot2")
library("reshape2")
library("varistran")
require(plyr)
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
    
    numericInput("filter", "", 10,),
    
    checkboxInput("varis", label = "Varistran Tansform Counts", value = T),
    
    checkboxInput("combine", label = "Combine by Replicates", value = F),
    
    textInput("file_name", label =  "Name of file to download", "file"),    
    
    uiOutput("choose_samples"),
    uiOutput("select_group"),
    downloadButton("downloadPlot", label = "Download Plot")
    
    
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
    get_plot_cols_c_v_t (count_info_table(), input$varis)
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
                         choices = 1:length(input$select_samples))
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
    if (input$combine == T){
      select_group_fun()      
    }
    
  })  
  output$choose_samples <- renderUI({
    checkboxGroupInput("select_samples", 
                       label = h5("Select which samples you would like to plot"),
                       choices = samples_from_counts(), 
                       selected = samples_from_counts()[1]) 
  })
  
  peak_tail_table <- reactive ({
    get_tail_vs_peak_pair(input$file_path)
  })
  
  
  output$data <- renderDataTable({
    if (input$select_plot_meth == 1){
      return(c_v_t_info_table())    
    }
    else if (input$select_plot_meth == 2){
      
      return(count_info_table())      
    }
    else{
      
      return(peak_tail_table()) 
      
    }
  })
  
  plot_calcs <- reactive({
    if (input$select_plot_meth == 1){
      plot <- make_plot_c_v_t(c_v_t_info_table(),input$select_samples, 
                              input$select_plot_meth, input$file_path, 
                              group_list(), input$combine)
      return(plot)    
    }
    else if (input$select_plot_meth == 2){
      return(count_info_table(c_v_t_info_table()))      
    }
    else{
      plot <- make_plot_c_v_t(peak_tail_table(), input$select_samples, 
                              input$select_plot_meth, input$file_path, 
                              group_list(), input$combine)
      
    }    
  }) 
  output$plot <- renderPlot({  
    plot_calcs()
  })
  output$downloadPlot <- downloadHandler(
    filename = function(){
      paste0(input$file_name, '.eps')
    },
    content = function(file){
      setEPS(width = 10)
      postscript(file)
      print(plot_calcs())    
      dev.off()       
    })  
  
}

make_info_frame <- function (tt_location){
  counts_csv <- read.grouped.table(paste0(tt_location, "/counts.csv"))  
  return (data.frame(counts_csv))
}

get_plot_cols_c_v_t <- function (count_table, varis){
  
  counts <- count_table [,grep("Count.*",colnames(count_table))]
  # Varistran or log2 counts depending on input
  if (varis == T){
    counts <- varistran::vst(counts)
  }
  else{
    counts <- log2(counts)
  }
  
  mean_tails <- count_table [,grep("Tail.*",colnames(count_table))]
  final <- data.frame(count_table$Annotation.gene, counts, mean_tails)
  return(final)
  
}

get_tail_vs_peak_pair <- function (file_path){
  peak_pair <- read.grouped.table(paste0(file_path, "/individual-pairs.csv"))  
  return (data.frame(peak_pair))
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


process_pp_t_frame <- function (cols_i_want, samples, 
                                title, group_list, combine){  
  for (sample in samples){
    peak_1 <- cols_i_want[grep(paste0("*", sample, "*peak1*"), colnames(cols_i_want))] 
    peak_2 <-cols_i_want[,grep(paste0("*", sample, "*peak2*"), colnames(cols_i_want))] 
  }
  samps <<- samples
  tit <<- title
  gloo <<- group_list
  
}

make_plot_c_v_t <- function (df, samples, name, title, group_list, combine){
  first <<- df
  if (name ==1){
    t <- paste ("Counts vs Tail Length",title)
    x <- "Count"
    y <- "Tail Length"
    multi <- ""
    names_col <- "count_table.Annotation.gene"
  }
  else if (name ==2){
    t <- paste ("Counts vs Peak Shift",title)
    x <- "Count"
    y <- "Peak Shift"
    multi <- ""
  }
  else{
    t <- paste ("Peak Shift vs Tail Length",title)
    x <- "Peak Shift"
    y <- "Tail Length Change"
    multi <- "*"
    names_col <- "Annotation.gene"
    
  }
  
  count_frame <- data.frame()
  count <- 1
  for (sample in samples){
    sample_col_count <- df[,grep(paste0("Count.", sample,multi), colnames(df))] 
    sample_col_tail <- df[,grep(paste0("Tail.", sample, multi), colnames(df))] 
    
    if (combine ==T){
      to_bind <- data.frame(df[,names_col], sample_col_count ,sample_col_tail, 
                            sample, group_list[[count]][1])      
    }
    else{
      to_bind <- data.frame(df[,names_col], sample_col_count ,sample_col_tail, 
                            sample)
    }
    count_frame <- rbind(count_frame,to_bind)
    count <- count+1
    
  }
  
  
  
  processed_frame_pp_v_t <- process_pp_t_frame (
    count_frame, samples,
    t, group_list, combine
  ) 
  if (combine == T){
    count_frame <- combine_by_reps(count_frame)   
    return(
      ggplot(data=count_frame, aes (x=mean_c, y=mean_t))+
        geom_point(aes(colour= group, group= group))+
        labs(title = t, x = x, y = y)
    )
  }
  else{    
    return(
      ggplot(data=count_frame, aes (x=sample_col_count, y=sample_col_tail))+
        geom_point(aes(colour= sample, group= sample))+
        labs(title = t, x = x, y = y)
    )
  }
}

combine_by_reps <- function(count_frame) {
  colnames(count_frame) <- c("gene", "count", "tail", "sample", "group")
  
  cdata <- ddply(count_frame, c("group", "gene"), summarise,
                 
                 mean_t = mean(tail,na.rm = T),
                 mean_c = mean(count, na.rm = T)
  )
  return(cdata)
}



shinyApp (ui = ui, server = server)
