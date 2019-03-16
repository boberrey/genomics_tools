# standard qc bargraph for ATAC results

library(ggplot2)
library(reshape2)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# table of values
read_stats <- read.delim(args[1])

# x column
x_name <- args[2]

# y column
y_name <- args[3]

# pdf name
out_name <- args[4]

qc_color = "#3A5795"

pdf(out_name,width=6,height=2.25 + 0.15*nrow(read_stats))

p1 <- (ggplot(read_stats, aes_string(x = x_name, y = y_name)) 
  + theme_bw() 
  + theme(panel.grid.major=element_blank(), panel.grid.minor= element_blank()) 
  + geom_bar(stat="identity",fill=qc_color) 
  + coord_flip()
)
p1
dev.off()