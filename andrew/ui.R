shinyUI(fluidPage(
  headerPanel('SCP2.0'),
  
  sidebarPanel(
    textInput("select_genes", label = h5("Select Your 
                                         Favourite Gene(s) or a Peak"),
              value = "Peak1"),
    checkboxInput("compare_datasets", "Compare Datasets?"),
    selectInput(inputId = 'dataset','Select Your Dataset', list_datasets()),
    conditionalPanel(
      condition = "input.compare_datasets == true",
      selectInput("data_set", "Data Set 2",
                  list_datasets())
    ),
    selectInput("gff_select", label = h4("GFF File Selection"), 
                choices = find_gff_files(gff_file_path), 
                selected =  find_gff_files(gff_file_path)[1]),
    checkboxGroupInput("select_bam_files", label = h4("Select the Relevant 
                                                      Bam FIles"),
                       choices = find_bam_files(bam_file_path), 
                       selected = find_bam_files(bam_file_path)[1]),br(),
    checkboxInput("merge", label = "Merge Replicates", value = F),
    conditionalPanel(
      condition = "input.merge == true",
      textInput("n_replicates", label = h5("How many replicates do you have?
                                           Enter an integer or format: 
                                           1:2 3:5 etc if replicate numbers
                                           are odd"),
                   value = 2)
    ),
    checkboxInput("legend", label = "Display Legend", value = F), br(),
    checkboxInput("show_reads", label = "Show Reads", value = T), br(),
    sliderInput("xslider", label= 'x axis slider', min=0, max=400,
                value =150, step = 25,ticks = TRUE, 
                sep = ","),
    sliderInput("ad_slider", label= 'number of adapter bases', min=0, max=23,
                value =0, step = 1,ticks = TRUE, 
                sep = ","),
    sliderInput("al_length", label= 'alignment length range', min=0, max=300,
                value =c(0,300)),
  h3("save"),
  downloadButton("downloadPlot", label = "Download Plot")  
),
  mainPanel(
    h2("Plot Output"),
    plotOutput('scp_plot'),
    textOutput('selected_dataset'),
    textOutput('read_files'),
    dataTableOutput("gff_rows"),
    conditionalPanel(
      condition = "input.show_reads == true",
      dataTableOutput('print_poly_a_counts')
    )
  )
))
