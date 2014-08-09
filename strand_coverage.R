# Plot strand coverage
# 3/13/2014 - D Sims
###############################################################################################################################
library("ggplot2")
library("reshape2")

Sys.time()
getwd()

args <- commandArgs( trailingOnly = TRUE )
if( length(args) != 2 ) stop( "ERROR: not enough args supplied" )

input_file <- args[1]
outdir <- args[2]

input_data <- data.frame( read.table( input_file, sep="\t", header=TRUE ) )
input_data$negReverse <- -input_data$Reverse
input_data$fwd_proportion<- input_data$Forward/(input_data$Forward+input_data$Reverse)
input_data<- within( input_data, Amplicon <- factor( Amplicon, levels=Amplicon[order(Median, decreasing=FALSE)]))

# plot all strand data
all.melt <- melt(input_data[, c("Amplicon", "Forward", "negReverse")], id.vars="Amplicon", variable.name="Direction", value.name="Coverage")

all_coverage_plot<- ggplot(all.melt, aes(Amplicon)) +
        geom_bar(data=all.melt[which(all.melt$Direction == "Forward"), ], aes(y=Coverage, fill="royalblue3"), stat="identity", position="identity") +
        geom_bar(data=all.melt[which(all.melt$Direction == "negReverse"), ], aes(y=Coverage, fill="red3"), stat="identity", position="identity") +
        theme_bw() +
        scale_fill_manual(values=c("red3", "royalblue3"), name="Direction", labels=c("Reverse", "Forward")) +
        ggtitle("Median Strand Coverage Distribution for All Amplicons") +
        theme(plot.title = element_text(size=16, face="bold")) +
        theme(axis.ticks.x = element_blank()) + 
        theme(axis.text.x = element_blank()) +
        theme(axis.text.y=element_text(size=12), axis.title=element_text(size=14, face="bold")) + 
        theme(legend.title=element_text(size=14, face="bold")) +
        theme(legend.text=element_text(size=12, face="bold"))


# Get a subset with a large variance for plotting below
high_variance <- subset( input_data, Forward/Reverse>1.95 | Reverse/Forward>1.95 )
high_variance <- within( high_variance, Amplicon <- factor( Amplicon, levels=Amplicon[order(fwd_proportion, decreasing=FALSE)]))
var.melt <- melt(high_variance[, c("Amplicon", "Forward", "negReverse")], id.vars="Amplicon", variable.name="Direction", value.name="Coverage")

high_variance_plot <- ggplot(var.melt, aes(Amplicon)) +
        geom_bar(data=var.melt[which(var.melt$Direction == "Forward"), ], aes(y=Coverage, fill="royalblue3"), stat="identity", position="identity") +
        geom_bar(data=var.melt[which(var.melt$Direction == "negReverse"), ], aes(y=Coverage, fill="red3"), stat="identity", position="identity") +
        theme_bw() +
        scale_fill_manual(values=c("red3", "royalblue3"), name="Direction", labels=c("Reverse", "Forward")) +
        ggtitle("Forward and Reverse Strand Proportion for High Variance Amplicons") +
        theme(plot.title = element_text(size=16, face="bold")) +
        theme(axis.text=element_text(size=12), axis.title=element_text(size=14, face="bold")) + 
        theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5, )) + 
        theme(legend.title=element_text(size=14, face="bold")) +
        theme(legend.text=element_text(size=12, face="bold"))

ggsave( filename=paste0( outdir, "/all_coverage_bias.png"), plot=all_coverage_plot, width=14, height=7)
ggsave( filename=paste0( outdir, "/strand_bias.png"), plot=high_variance_plot, width=14, height=7)
