
shinyUI(fluidPage(  
  # Hide all error messages
  tags$style(type="text/css",
             ".shiny-output-error { visibility: hidden; }",
             ".shiny-output-error:before { visibility: hidden; }"
  ),
  titlePanel("Poly(A) Plotter"),
  em(helpText("created by Andrew Pattison, Jack Xu and Paul Harrison for the Beilharz Lab", align = "right")),
  helpText("
           This app analyses the percentage population against poly-A tail length in a 
           specific gene or peak, from the selected sample files.It also provides a range
           of other potentially useful information about a PAT-Seq or Re-PAT"), br(),
  
  sidebarPanel(
    
    uiOutput("select_file_path"),
    textInput("select_genes", label = ("Select your 
                                         favourite gene(s) or a peak(s) 
                                         separated by a space"),
              value = "Peak1873"),
    uiOutput("gff_files"),
    checkboxInput("alt_plot", label = "View information on the sequences preceeding the 
                  poly (A)-tail", value = F),
    conditionalPanel(
      condition = "input.alt_plot == true",
      checkboxInput("poly_a_pileup", label = "Show raw read and poly (A)-tail lengths", value = F),
      conditionalPanel(
        condition = "input.poly_a_pileup == true",
        checkboxInput("order_alt", label = "Order reads by alignment length", value = T)  
      )
    ),
    conditionalPanel(
      condition = "input.alt_plot == true && input.poly_a_pileup == false",      
      checkboxInput("alt_cumu_dis", label = "Set maximum reads to 100%", value = T)
      
    ),
    checkboxInput("all_reads", label = "Include reads that do not have a 
                  poly (A)-tail", value = F),
    checkboxInput("merge", label = "Combine samples", value = F),
    conditionalPanel(
      condition = "input.merge == false",
      uiOutput("bam_files") 
      
    ),
    
    conditionalPanel(
      condition = "input.merge == true",   
      uiOutput("select_group")
    ),
    ("save"),
    downloadButton("downloadPlot", label = "Download plot"),
    
    checkboxInput("legend", label = "Display a legend on the plot", value = T),
    checkboxInput("mm_frame", label = "Show poly (A)-tail length means and medians from these reads", value = T),  
    checkboxInput("gff_info", label = "Display info about the genomic regions used in this plot", value = T),
    checkboxInput("show_reads", label = "Show the reads that were used in making this plot", value = T),     
    sliderInput("xslider", label= 'x axis slider', min=0, max=400,
                value =c(0, 300), step = 25,ticks = TRUE, 
                sep = ","),
    sliderInput("ad_slider", label= 'number of adapter bases (filters out reads that
                were not sequenced completely to the end of the poly (A)-tail', min=0, max=23,
                value =0, step = 1,ticks = TRUE, 
                sep = ","),
    sliderInput("al_length", label= 'alignment length range (gets reads that had a 
                genomic aligment within the selected length)', min=0, max=400,
                value =c(0,400))
    
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
