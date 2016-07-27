# Create a length vs median coverage scatterplot for AmpliconCoverageAnalysis Plugin
#
# 3/3/2014 - D Sims
######################################################################################################################
getwd()
#suppressPackageStartupMessages()

# Get some commandline args
args <- commandArgs( trailingOnly = TRUE )
if( length( args ) != 3 ) stop( "Script requires 3 arguments: sample_name, minimum_coverage_threshold, and output dir." )
sample    <- args[1]
threshold <- args[2]
outdir    <- args[3]

# Make sure we have ggplot2 installed
if(!require(ggplot2)) stop( "Need to install ggplot2" )

coverage_data <- data.frame(read.table(file=paste0(outdir,"/AllAmpliconsCoverage.tsv"), sep="\t", header = TRUE))

scatter <- ggplot(coverage_data, aes(x = Length, y = Median)) +
           geom_point() +
           geom_hline( yintercept = as.numeric(threshold), linetype = "dashed", color = "red3" ) +
           geom_text(aes(0, as.numeric(threshold), label=paste0(threshold, "X"), vjust = -1), size=4, fontface="bold", color="dodgerblue4") + 
           theme_bw() +
           theme(
               panel.background = element_rect( fill = "aliceblue"),
               panel.grid.major = element_line( color = "lightcyan4", size=0.2, linetype="dotted" ),
               panel.grid.minor = element_line( color = "lightcyan4", size=0.2, linetype="dotted" ),
               axis.title = element_text( size=14, face="bold" ),
               plot.title = element_text( size=16, face="bold")
           ) + 
           xlab( "Length of Amplicon" ) +
           ylab( "Median Coverage for Amplicon" ) +
           ggtitle( paste0( "Coverage Versus Amplicon Length Plot for\n", sample ) )

ggsave( filename = paste0( outdir, "/Amp_Coverage_vs_Length_Plot.png" ), plot = scatter )
ggsave( filename = paste0( outdir, "/Amp_Coverage_vs_Length_Plot.pdf" ), plot = scatter )
