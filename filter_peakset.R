# Filter peakset

suppressPackageStartupMessages({
  library(GenomicRanges)
})

# Read in command line arguments
args <- commandArgs(trailingOnly = TRUE)

peak_df <- read.delim(args[1], header = FALSE)

names(peak_df) <- c("chrom", "start", "end", "peak_name", "macs2_score", "spm_score")

# Make GenomicRanges for merged peaks
peak_ranges <- with(peak_df, GRanges(seqnames = chrom, 
                                     IRanges(start, end), 
                                     score1 = macs2_score, 
                                     score2 = spm_score,
                                     id = peak_name))

# Keep only standard chromosomes
standard_peaks <- keepStandardChromosomes(peak_ranges, pruning.mode = "tidy")

# Remove chromosomeY (because everyone else is doing it?)
no_Y <- standard_peaks[seqnames(standard_peaks) != "chrY"]

# Save filtered peak file
df_to_save <- data.frame(no_Y)
df_to_save <- df_to_save[,c("seqnames", "start", "end", "id", "score1", "score2")]

# Write raw counts matrix
write.table(df_to_save, file = args[2], 
            sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)

