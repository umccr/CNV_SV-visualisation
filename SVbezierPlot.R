################################################################################
#
#   Description: Script generating a bezier curves-like plot representing structural variants listed in specified BEDPE file. This script catches the arguments from the command line and pass them to the SVbezierPlot.Rmd script to produce the report and generate a bezier curves-like plot representing structural variants.
#   NOTE: currently the script plots take into account only the start of each segement (second and fifth columns in the BEDPE file) to plot the bezier curves between the different genomic locations of individual events.
#
#   Command line use example: Rscript SVbezierPlot.R --bedpe example_data/structural/example-manta-pass.bedpe --output example_data/structural/SVbezierPlot
#
#   bedpe:     BEDPE file to be processed
#   output:    The output plot name. Note, the plot will be saved as png file. If no output plot name is specified the file will have the same name as the input bedpe file with added "_SVbezier" at the end and the extension changes to ".png"
#
################################################################################

#===============================================================================
#    Catching the arguments
#===============================================================================

suppressMessages(library(optparse))
suppressMessages(library(rmarkdown))

option_list <- list(
  make_option(c("-i", "--bedpe"), action="store", default=NA, type='character',
              help="BEDPE file to be processed"),
  make_option(c("-o", "--output"), action="store", default=NA, type='character',
              help="Output plot name")
)

opt <- parse_args(OptionParser(option_list=option_list))

bedpe.file <- opt$bedpe
output.plot <- opt$output

##### Change the file extension to .maf if output file name is not specified
if ( is.na(output.plot) ) {
  	output.plot <- unlist(strsplit(bedpe.file, split='.', fixed=TRUE))
  	output.plot <- paste0(paste(output.plot[-length(output.plot)], collapse = '.'), "_SVbezier.png")

} else {
    output.plot <- paste0(output.plot, ".png")
}

##### Run the
rmarkdown::render(input="SVbezierPlot.Rmd", params = list(bedpe = bedpe.file, outfile = output.plot))
