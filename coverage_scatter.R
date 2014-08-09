# Create a length vs median coverage scatterplot for AmpliconCoverageAnalysis Plugin
#
# 3/3/2014 - D Sims
######################################################################################################################

Sys.time()
getwd()

# Get some commandline args
args <- commandArgs( trailingOnly = TRUE )
if( length( args ) !=3 ) stop( "Script requires 3 arguments: sample_name, minimum_coverage_threshold, and output dir." )
sample <- args[1]
threshold <- args[2]
outdir <- args[3]

# Make sure we have ggplot2 installed
if(!require(ggplot2)) stop( "Need to install ggplot2" )

file <- paste0( outdir, "/AllAmpliconsCoverage.tsv" )
coverage_data <- read.table( file = file, sep = "\t", header = TRUE ) 

scatter <- ggplot( coverage_data, aes( x = Length, y = Median ) ) +
           theme_bw() +
           #theme(panel.background = element_rect( fill = "lightcyan2")) +
           theme(panel.background = element_rect( fill = "aliceblue")) +
           #theme(panel.background = element_rect( fill = "azure2")) +
           theme(panel.grid.major = element_line( colour = "lightcyan4", size=0.2, linetype="dotted" )) +
           theme(panel.grid.minor = element_line( colour = "lightcyan4", size=0.2, linetype="dotted" )) +
           geom_point() +
           geom_hline( yintercept = as.numeric(threshold), linetype = "dashed", color = "red3" ) +
           geom_text(aes( 0, as.numeric(threshold), label=paste0(threshold, "X"), vjust = -1), size=4, face="bold", color="dodgerblue4") +
           theme(axis.title = element_text( size=14, face="bold" )) +
           xlab( "Length of Amplicon" ) +
           ylab( "Median Coverage for Amplicon" ) +
           theme(plot.title = element_text( size=16, face="bold")) +
           ggtitle( paste0( "Coverage Versus Amplicon Length Plot for\n", sample ) )

ggsave( filename = paste0( outdir, "/Amp_Coverage_vs_Length_Plot.png" ), plot = scatter )
ggsave( filename = paste0( outdir, "/Amp_Coverage_vs_Length_Plot.pdf" ), plot = scatter )
