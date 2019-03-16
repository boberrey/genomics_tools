# Plot the insert size distribution of a paired-end sequencing library

# Usage:
# get_insert_size_dist.R <bam_file> <plot_filename> <window_limit>

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomicAlignments)
  library(ggplot2)
})

# Read in command line arguments
args <- commandArgs(trailingOnly = TRUE)
bam_file <- BamFile(args[1])
output_file <- args[2]
window_size <- as.numeric(args[3])


# Function to read bam file and get fragment lengths
bam_to_length_table <- function(bam_file){
  reads <- as(readGAlignmentPairs(bam_file), "GRanges")
  tbl <- table(width(reads))
  df <- as.data.frame(tbl)
  names(df) <- c("insert_size", "fraction")
  # Convert position from factor
  df$insert_size <- as.numeric(as.character(df$insert_size))
  df$fraction <- df$fraction/sum(df$fraction)
  return(df)
}

data_df <- bam_to_length_table(bam_file)

p <- (
  ggplot(data = data_df, aes(x=insert_size, y=fraction))
  + geom_line(color="red", size=0.8)
  + scale_x_continuous(limits = c(0, window_size), 
                       breaks = seq(0, window_size, window_size/4), 
                       expand = c(0,0))
  + scale_y_continuous(expand = c(0.0001,0.0001, 0.0005, 0.0005))
  + xlab("Insertion size")
  + ylab("Fraction inserts")
  + theme_bw()
  + theme(panel.grid.major=element_blank(), 
          panel.grid.minor= element_blank(), 
          plot.margin = unit(c(0.25,1,0.25,1), "cm")) 
  
)

aspect_ratio <- 1.5
height <- 3
ggsave(output_file, height = height, width = aspect_ratio*height)

