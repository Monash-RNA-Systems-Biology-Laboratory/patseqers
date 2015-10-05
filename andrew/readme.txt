This is a Shiny App. It can be used to plot an inverse cumulative distribution of BAM files from a Tail-Tools output.


For info on Tail-Tools see:
https://github.com/Victorian-Bioinformatics-Consortium/tail-tools

Set up:

Requires Rstudio to be installed.

Requires the following R packages to be installed:

Shiny
Rsamtools

To install Rsamtools simply type into R.

source("http://bioconductor.org/biocLite.R")
biocLite("Rsamtools")

For info on RStudio/installation see:
https://www.rstudio.com/

For info on Shiny/installation see:
http://shiny.rstudio.com/

To make the download plot button work:
extrafont needs to be installed and set up.

install.packages("extrafont")
This command only works in terminal R for some reason.
library("extrafont")
font_import()
loadfonts(device = "postscript")

To use:
Set up shiny server, or run locally from RStudio.

Place or symlink (ln -s path/to/dir) to directories containing relevant bam and gff files into the directory containing Server.R, ui.R and major_calclations.R. You can then select datasets from the app based on directories.
