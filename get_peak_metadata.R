# Get peak information to use later on

library(RMariaDB)
library(GenomicFeatures)
library(GenomicRanges)
library(org.Mm.eg.db)

mm10.refseq.db <- makeTxDbFromUCSC(genome = "mm10", tablename = "ncbiRefSeqCurated")

peak_df <- read.delim('output/peaks/merged_peaks/filtered_peak_set.bed', header = FALSE)

names(peak_df) <- c("chrom", "start", "end", "peak_name", "macs2_score", "spm_score")


# Just use the nice built-in's...
chrs <- peak_df$chrom
locus <- peak_df$start + (peak_df$end - peak_df$start)/2
nearest_tss <- nearestTSS(chrs, locus, species = "Mm")

updated_df <- cbind(peak_df, nearest_tss)

