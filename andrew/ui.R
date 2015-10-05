
shinyUI(fluidPage(  
  # Hide all error messages
#   tags$style(type="text/css",
#              ".shiny-output-error { visibility: hidden; }",
#              ".shiny-output-error:before { visibility: hidden; }"
#   ),
  titlePanel("Poly(A) Plotter"),
  em(helpText("created by Andrew Pattison, Jack Xu and Paul Harrison for the Beilharz Lab", align = "right")),
  helpText("
           This app analyses the percentage population against poly-A tail length in a 
           specific gene or peak, from the selected sample files."), br(),
  
  sidebarPanel(
    
    uiOutput("select_file_path"),
    textInput("select_genes", label = h5("Select Your 
                                         Favourite Gene(s) or a Peak(s) 
                                         Separated by a Space"),
              value = "Peak1873"),
    uiOutput("gff_files"),
    checkboxInput("alt_plot", label = "Plot Genomic Alignments", value = F),
    conditionalPanel(
      condition = "input.alt_plot == true",
      checkboxInput("poly_a_pileup", label = "Pile up Reads", value = F),
      checkboxInput("order_alt", label = "Order Reads By Width and then Poly (A) Tail Length", value = T)       
    ),
      conditionalPanel(
        condition = "input.poly_a_pileup == false",      
      checkboxInput("alt_cumu_dis", label = "Plot a Cumulative Distribution of Read Lengths", value = T)
      
    ),
    checkboxInput("all_reads", label = "Include Non-Poly (A) Reads", value = F),
    checkboxInput("merge", label = "Merge Replicates", value = F),
    conditionalPanel(
      condition = "input.merge == false",
      uiOutput("bam_files") 
      
    ),
    
    conditionalPanel(
      condition = "input.merge == true",   
      uiOutput("select_group")
    ),
    h3("save"),
    downloadButton("downloadPlot", label = "Download Plot"),
    
    checkboxInput("legend", label = "Display Legend", value = T),
    checkboxInput("mm_frame", label = "Show Poly-(A) Tail Length Means and Medians from these Reads", value = T),  
    checkboxInput("gff_info", label = "Display GFF Info", value = T),
    checkboxInput("show_reads", label = "Show Reads", value = T),     
    sliderInput("xslider", label= 'x axis slider', min=0, max=400,
                value =c(0, 300), step = 25,ticks = TRUE, 
                sep = ","),
    sliderInput("ad_slider", label= 'number of adapter bases', min=0, max=23,
                value =0, step = 1,ticks = TRUE, 
                sep = ","),
    sliderInput("al_length", label= 'alignment length range', min=0, max=300,
                value =c(0,300))
    
  ),
  mainPanel(
    tabsetPanel(
      tabPanel("Plot",
               
               h2("Plot Output",  height = "600"),   
               plotOutput('scp_plot', brush = "plot_brush"),
               textOutput('frame_type'),
               dataTableOutput("seq_plot"),
               em(textOutput("gene_info")),
               textOutput('read_files'),
               textOutput('groups'),
               textOutput('n_reps')
      ),
      tabPanel("Info",
               conditionalPanel(
                 condition = "input.gff_info == true",
                 dataTableOutput("gff_rows")
               ),
               conditionalPanel(
                 condition = "input.mm_frame == true",
                 dataTableOutput("means_frame")
               ),
               conditionalPanel(
                 condition = "input.show_reads == true",
                 dataTableOutput('print_poly_a_counts')
               )
      )
    )
  )
))
