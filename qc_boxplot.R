# standard qc bargraph for ATAC results

library(ggplot2)
library(reshape2)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# output from pipeline
read_stats <- read.delim(args[1])

# axis label
label <- args[2]

# pdf name
out_name <- args[3]

# Melt table and assign column headers to be 'metrics'
melted_stats <- melt(read_stats, var="metric")

# convert metrics to factors
melted_stats$metric <- factor(melted_stats$metric, levels=unique(melted_stats$metric))

qc_color = "#3A5795"

pdf(out_name)
(ggplot(melted_stats, aes(metric, value))
  + theme_bw()
  + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
  + geom_boxplot(color=qc_color) 
  + coord_flip() 
  + ylab(label)
)
dev.off()
  
