This is a Shiny App. It can be used to plot an inverse cumulative distribution of BAM files from a Tail-Tools output.


For info on Tail-Tools see:
https://github.com/Victorian-Bioinformatics-Consortium/tail-tools

Set up:

Requires the following R packages to be installed:

jsonlite
Rstudio
Shiny
Rsamtools 

To install Rsamtools simply type into R.

To install jsonlite
install.packages('jsonlite')

source("http://bioconductor.org/biocLite.R")
biocLite("Rsamtools")

For info on RStudio/installation see:
https://www.rstudio.com/

For info on Shiny/installation see:
http://shiny.rstudio.com/

To use:
Set up shiny server, or run locally from RStudio. 

Place directories containing relevant bam and gff files into the directory containing Server.R, ui.R and major_calclations.R. You can then select datasets from the app based on directories. Additionally, a symbolic link (ln -s <path to dir >) can be made to a tail-tools output from the directory containing the app and the application will be able to find the files it requires. 
