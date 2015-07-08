shinyUI(fluidPage(
  titlePanel("Poly(A) Plotter"),
  em(helpText("created by Andrew Pattison, Jack Xu and Paul Harrison for the Beilharz Lab", align = "right")),
  helpText("
           This app analyses the percentage population against poly-A tail length in a 
           specific gene or peak, from the selected sample files."), br(),
  
  sidebarPanel(
    textInput("select_genes", label = h5("Select Your 
                                         Favourite Gene(s) or a Peak"),
              value = "Peak1"),
    selectInput(inputId = 'dataset','Select Your Dataset', list_datasets()),
    checkboxInput("compare_datasets", "Compare Datasets"),
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
    h3("save"),
    downloadButton("downloadPlot", label = "Download Plot"),
    checkboxInput("merge", label = "Merge Replicates", value = F),
    conditionalPanel(
      condition = "input.merge == true",
      textInput("n_replicates", label = h5("How many replicates do you have?
                                           Enter an integer or format: 
                                           1:2 3:5 etc if replicate numbers
                                           are odd"),
                   value = 2)
    ),
    checkboxInput("legend", label = "Display Legend", value = T),
    checkboxInput("gff_info", label = "Display GFF Info", value = F),
    checkboxInput("show_reads", label = "Show Reads", value = F),     
    sliderInput("xslider", label= 'x axis slider', min=0, max=400,
                value =c(0, 150), step = 25,ticks = TRUE, 
                sep = ","),
    sliderInput("ad_slider", label= 'number of adapter bases', min=0, max=23,
                value =0, step = 1,ticks = TRUE, 
                sep = ","),
    sliderInput("al_length", label= 'alignment length range', min=0, max=300,
                value =c(0,300))

),
  mainPanel(
    h2("Plot Output",  height = "600"),
    plotOutput('scp_plot'),
    em(textOutput("gene_info")),
    textOutput('read_files'),
    conditionalPanel(
    condition = "input.gff_info == true",
    dataTableOutput("gff_rows")
    ),
    conditionalPanel(
      condition = "input.show_reads == true",
      dataTableOutput('print_poly_a_counts')
    )
  )
))
