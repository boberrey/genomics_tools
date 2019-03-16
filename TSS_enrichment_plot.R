# Make TSS enrichment plots

# Usage:
# TSS_enrichment_plot.R <bam_file_path> <output_file_path> <genome> <window_size> <n_cores>

# Ben Ober-Reynolds

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomicAlignments)
  library(RMariaDB)
  library(GenomicFeatures)
  library(preprocessCore)
  library(ggplot2)
})


# Read in command line arguments
args <- commandArgs(trailingOnly = TRUE)
bam_filepath <- args[1]
output_filepath <- args[2]
gnome <- args[3]
window_size <- as.numeric(args[4])
n_cores <- as.numeric(args[5])

# # Read in command line arguments (debug)
# bam_filepath <- "output/bams/deduped/"
# output_filepath <- "output/plots/TSS/"
# gnome <- "mm10"
# window_size <- 2000
# n_cores <- 10


# Gather TSS's
if(gnome == "mm10"){
  gnome.refseq.db = makeTxDbFromUCSC(genome = "mm10", tablename = "ncbiRefSeqCurated")
} else {
  gnome.refseq.db = makeTxDbFromUCSC(genome = "hg38", tablename = "ncbiRefSeqCurated")
}

# Get non-redundant TSS sites
tss <- unique(
  keepStandardChromosomes(
    resize(transcripts(gnome.refseq.db), width = 1, fix = 'start'), pruning.mode = "tidy"))

# Remove mito
tss <- tss[seqnames(tss) != "chrM"]

# Function to read bam file into Tn5-shifted genomicranges
bam_to_shifted_range <- function(bam_file){
  reads <- as(readGAlignments(bam_file), "GRanges")
  # Resize to isolate the ends of fragments (the 5' end of each read)
  reads <- resize(reads, 1, fix = 'start')
  # Adjust insertion based on strand
  # ('+' stranded +4 bp, '-' stranded -5 bp)
  reads <- shift(reads, ifelse(strand(reads) == '+', 4, -5))
  reads
}

# Moving average function:
moving_average = function(x, n = smooth_size){filter(x, rep(1/n, n), sides=2)}

# Insertion window function
get_insertion_window <- function(bam_file, query_range, window_size, smooth_window=51){
  # Read in bam file and find overlaps
  read_range <- bam_to_shifted_range(bam_file)
  overlaps <- findOverlapPairs(read_range, query_range, 
                               ignore.strand = TRUE, maxgap = window_size-1)
  # Correct for TSS directionality
  flips <- ifelse(strand(second(overlaps)) == "-", -1, 1)
  diffs <- (start(first(overlaps)) - start(second(overlaps)))*flips
  dist_df <- as.data.frame(table(diffs))
  # Convert position from factor
  dist_df$diffs <- as.numeric(as.character(dist_df$diffs))
  # Normalize to mean of counts from positions +/- 1900-2000
  norm_factor <- mean(c(head(dist_df$Freq, 100), tail(dist_df$Freq, 100)))
  dist_df$Freq <- dist_df$Freq / norm_factor
  # Get smoothing line
  dist_df$smoothed <- moving_average(dist_df$Freq, smooth_window)
  dist_df
}

# Function to plot TSS enrichment
plot_and_save_TSS <- function(df, output_filepath, sample_name, window_size){
  p <- (
    ggplot(df, aes(x=diffs, y=Freq))
    + theme_bw()
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm")) 
    + geom_point(color="grey", size=1)
    + geom_line(aes(x=diffs, y=smoothed), color="red") 
    + scale_x_continuous(limits = c(-window_size, window_size), 
                         breaks = seq(-window_size, window_size, window_size/4), 
                         expand = c(0,0))
    + xlab("Position Relative to TSS")
    + ylab("Insertions")
  )
  aspect_ratio <- 1.5
  height <- 5
  ggsave(paste(output_filepath, "/", sample_name, "_TSS_enrichment.pdf", sep = ""), 
         height = height, width = aspect_ratio*height)
}


# Locate bam files
bam_files <- list.files(path = bam_filepath, pattern = "\\.bam$", full.names = TRUE)

# Prepare for use with bioconductor I/O tools
sample_names <- sapply(strsplit(bam_files, '\\.'), function(l) l[[1]])
sample_names <- basename(sample_names)
bam_file_list <- lapply(bam_files, BamFile)

# get all insertion window data frames
window_dfs <- mclapply(bam_files, function(x) get_insertion_window(x, tss, window_size), mc.cores = n_cores)

# Generate all plots
for(i in 1:length(window_dfs)){
  plot_and_save_TSS(window_dfs[[i]], output_filepath, sample_names[i], window_size)
}


# Write TSS enrichments to file
get_max_enrichment <- function(df, tss_window=50){
  round(max(df[(df$diffs > -tss_window) & (df$diffs < tss_window),]$smoothed),5)
}
enrichments <- lapply(window_dfs, get_max_enrichment)
enrich_df <- data.frame(as.matrix(enrichments))
enrich_df <- apply(enrich_df, 2, as.character)
rownames(enrich_df) <- sample_names
write.table(enrich_df, file = paste(output_filepath, "/", "TSS_enrichments.txt", sep=""), 
                   quote=FALSE, sep="\t", col.names = FALSE, row.names = TRUE)

