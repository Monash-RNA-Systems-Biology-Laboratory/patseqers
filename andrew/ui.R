shinyUI(fluidPage(
  titlePanel("Poly(A) Plotter"),
  em(helpText("created by Andrew Pattison, Jack Xu and Paul Harrison for the Beilharz Lab", align = "right")),
  helpText("
           This app analyses the percentage population against poly-A tail length in a 
           specific gene or peak, from the selected sample files."), br(),
  
  sidebarPanel(
    selectInput(inputId = 'dataset','Select Your Dataset', list_datasets()),
    checkboxInput("compare", label = "Compare Datasets", value = F),
    conditionalPanel(
      condition = "input.compare == true",
      selectInput(inputId = 'dataset2','Select Your Second Dataset', list_datasets())
    ),      
    textInput("select_genes", label = h5("Select Your 
                                         Favourite Gene(s) or a Peak(s) 
                                         Separated by a Space"),
              value = "Peak1"),
    selectInput("gff_select", label = h4("GFF File Selection"), 
                choices = find_gff_files(gff_file_path), 
                selected =  find_gff_files(gff_file_path)[1]),
    checkboxInput("merge", label = "Merge Replicates", value = F),
    conditionalPanel(
    condition = "input.merge == false",
    checkboxGroupInput("select_bam_files", label = h4("Select the Relevant 
                                                      Bam FIles"),
                       choices = find_bam_files(bam_file_path), 
                       selected = find_bam_files(bam_file_path)[1])
    ),
    
    conditionalPanel(
      condition = "input.merge == true",   
      uiOutput("select_group")
     # lapply(find_bam_files(bam_file_path), function(i) {
     # selectInput(paste0('snumber', i),  h5(paste0('Select a group for ', i)),
      #              choices = 1:5)
   #   })
      ),
    h3("save"),
    downloadButton("downloadPlot", label = "Download Plot"),
    
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
    textOutput('groups'),
    textOutput('n_reps'),
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
