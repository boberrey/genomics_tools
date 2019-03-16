# Create a normalized counts matrix for bulk ATAC data

# Ben Ober-Reynolds

library(GenomicRanges)
library(GenomicAlignments)
library(edgeR)
library(preprocessCore)
library(ggplot2)


#args <- commandArgs(trailingOnly = TRUE)

peak_df <- read.delim("output/peaks/merged/merged_peaks_fix_width.bed", header = FALSE)
#peak_df <- read.delim(args[1], header=FALSE)

names(peak_df) <- c("chrom", "start", "end", "peak_name", "macs2_score", "spm_score")

# Make GenomicRanges for merged peaks
peak_ranges <- with(peak_df, GRanges(seqnames = chrom, 
                                     IRanges(start, end), 
                                     score = spm_score, 
                                     id = peak_name))



# Locate bam files
#bam_file_path <- args[2]
bam_file_path <- "output/bams/deduped/"
bam_files <- list.files(path = bam_file_path, pattern = "\\.bam$", full.names = TRUE)

# Prepare for use with I/O tools
headers <- sapply(strsplit(bam_files, '\\.'), function(l) l[[1]])
headers <- sapply(strsplit(headers, "/"), function(l) l[[length(l)]])
bam_file_list <- lapply(bam_files, BamFile)
names(bam_file_list) <- headers


# Function to get counts in peaks
get_counts <- function(bam_file, peak_ranges){
  # Read in bam file as GenomicRanges
  reads <- as(readGAlignments(bam_file), "GRanges")
  # Resize to isolate the ends of fragments (the 5' end of each read)
  reads <- resize(reads, 1)
  # Adjust insertion based on strand
  # ('+' stranded +4 bp, '-' stranded -5 bp)
  reads <- shift(reads, ifelse(strand(reads) == '+', 4, -5))
  # Get overlaps
  countOverlaps(peak_ranges, reads, ignore.strand = TRUE)
}

# Get overlaps for each bam file
overlaps <- lapply(bam_file_list, function(x) get_counts(x, peak_ranges))

names(overlaps) <- headers
overlap_matrix <- do.call(cbind, overlaps)

# Normalize count matrix

# First get counts per million
cpm_tst <- cpm(overlap_matrix, log = TRUE, prior.count = 3)
quant_tst <- normalize.quantiles(cpm_tst)


# Get most variable regions
counts_df <- data.frame(quant_tst)
names(counts_df) <- headers
counts_df$mean_counts <- apply(quant_tst, 1, mean)
counts_df$sd_counts <- apply(quant_tst, 1, sd)

counts_df <- round(counts_df, digits = 2)

full_df <- cbind(peak_df, counts_df)
write.table(full_df, file = "tst_counts_table.txt", 
            sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
mean_5_16 <- apply(full_df[c("1p5pct25Hex1_S20", "1p5pct25Hex2_S25")], 1, mean)
mean_mES <- apply(full_df[c("mES_WT_1_S21", "mES_WT_2_S26")], 1, mean)
mean_mean <- apply(cbind(mean_5_16, mean_mES), 1, mean)
diff <- mean_5_16 - mean_mES


tst_df <- data.frame(cbind(mean_mean, diff))

pdf("tst_plot.pdf",width=4.2,height=3.3)

p <- qplot(x = c(1:10), y = c(1:10))
p
dev.off()
