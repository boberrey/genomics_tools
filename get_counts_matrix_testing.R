# Create a normalized counts matrix for bulk ATAC data

# Usage:
# get_counts_matrix.R <peak_bed_file> <bam_file_dir> <output_filename> <n_cores> <genome> <normalization>

# Ben Ober-Reynolds

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomicAlignments)
  library(preprocessCore)
  library(AnnotationDbi)
  library(edgeR)
})


# # Read in command line arguments
# args <- commandArgs(trailingOnly = TRUE)
# peak_file <- args[1]
# bam_file_path <- args[2]
# output_filename <- args[3]
# n_cores <- as.numeric(args[4])
# g <- args[5] # {'mm10', 'hg38'}
# normalization <- args[6] # {'peaks', 'housekeeping', 'none'}

# Read in command line arguments (DEBUG)
#args <- commandArgs(trailingOnly = TRUE)
setwd("/raid/USRdirs/ben/phase_sep/hexanediol_ATAC/20190207_GM_and_HL_deep/")
peak_file <- "output/peaks/processed/filtered_fix_width_peaks.bed"
bam_file_path <- "output/bams/deduped/"
output_filename <- "output/counts_matrix/norm_testing_count_matrix.txt"
n_cores <- as.numeric('10')
g <- "hg38"
normalization <- "housekeeping" # {'peaks', 'housekeeping', 'none'}


# Get correct genome
if(g == "mm10"){
  library(org.Mm.eg.db)
  sp <- "Mm"
} else if(g == "hg38"){
  library(org.Hs.eg.db)
  sp <- "Hs"
} else {
  print("Invalid genome. Quitting...")
  quit()
}

between <- function(x, low, high){return((x > low) & (x < high))}

# Get the merged and filtered peak set to get counts from
peak_df <- read.delim(peak_file, header=FALSE)

names(peak_df) <- c("chrom", "start", "end", "peak_name", "macs2_score", "spm_score")

# Make GenomicRanges for merged peaks
peak_ranges <- with(peak_df, GRanges(seqnames = chrom, 
                                     IRanges(start, end), 
                                     score = spm_score, 
                                     id = peak_name))

# Locate bam files
bam_files <- list.files(path = bam_file_path, pattern = "\\.bam$", full.names = TRUE)

# Prepare for use with bioconductor I/O tools
headers <- sapply(strsplit(bam_files, '\\.'), function(l) l[[1]])
headers <- basename(headers)
bam_file_list <- lapply(bam_files, BamFile)
names(bam_file_list) <- headers

# Specify parameters for reading in Bam files:
bam_param <- ScanBamParam(mapqFilter = 30)


# Function to get counts in peaks
get_counts <- function(bam_file, peak_ranges, bam_param){
  # Read in bam file as GenomicRanges
  reads <- as(readGAlignments(bam_file, param = bam_param), "GRanges")
  # Resize to isolate the ends of fragments (the 5' end of each read)
  reads <- resize(reads, 1, fix = 'start')
  # Adjust insertion based on strand
  # ('+' stranded +4 bp, '-' stranded -5 bp)
  reads <- shift(reads, ifelse(strand(reads) == '+', 4, -5))
  # Get overlaps
  countOverlaps(peak_ranges, reads, ignore.strand = TRUE)
}

# Get overlaps for each bam file (in parallel)
overlaps <- mclapply(bam_file_list, function(x) get_counts(x, peak_ranges, bam_param), mc.cores = n_cores)

names(overlaps) <- headers
overlap_matrix <- do.call(cbind, overlaps)


# Filter peaks that have 10x the 95th percentile?
#gmeans <- apply(overlap_matrix, 1, function(x){exp(mean(log(x + 1)))})


# Perform normalization

if(normalization == 'housekeeping'){
  ### Housekeeping Normalization ###
  # Normalize to the least variable 5000 promoter-associated peaks
  # In the TCGA paper they define 'promoter regions' as being within -1000 to +100
  # window around a TSS
  # The 'nearestTSS' function reports positive distances as being 'upstream' of the TSS
  # while negative distances are in the gene body. Which is the one other people use?
  nearest_tss <- nearestTSS(peak_df$chrom, 
                            peak_df$start + (peak_df$end - peak_df$start)/2, 
                            species = sp)
  
  tss_df <- cbind(peak_df[,1:4], "distance" = nearest_tss$distance, overlap_matrix)
  promoter_df <- tss_df[between(tss_df$distance, -1000, 100),]
  # Perform cpm+quantile normalization on this set 
  #normed_m <- normalize.quantiles(cpm(promoter_df[,6:ncol(promoter_df)], log=TRUE, prior.count=5))
  normed_m <- normalize.quantiles(cpm(promoter_df[,6:ncol(promoter_df)], log=TRUE, prior.count=5))
  normed_df <- as.data.frame(cbind(promoter_df[,1:5], normed_m))
  #normed_df$CV <- apply(normed_df[,6:ncol(normed_df)], 1, function(x){sd(x)/mean(x)})
  #norm_peaks <- head(normed_df[order(normed_df$CV, decreasing = F),], 4000)
  
  # Just force the median?
  normed_df$CV <- apply(normed_df[,6:ncol(promoter_df)], 1, function(x){sd(x)/mean(x)})
  norm_peaks <- normed_df[order(normed_df$CV, decreasing = T),]
  
  # Not pick out these peaks from the raw counts
  rownames(promoter_df) <- as.character(promoter_df$peak_name)
  norm_peaks <- promoter_df[as.character(norm_peaks$peak_name),]
  # # filter any peaks that have less than 5 reads in any sample
  # promoter_df <- promoter_df[apply(promoter_df[,6:ncol(promoter_df)], 1, function(x){all(x > 5)}),]
  # # How to define 'least variable' peaks? For now just use CV.
  # promoter_df$CV <- apply(promoter_df[,6:ncol(promoter_df)], 1, function(x){sd(x)/mean(x)})
  # norm_peaks <- head(promoter_df[order(promoter_df$CV, decreasing = F),], 4000)
  # Do a 'DEseq-like' median of ratios approach to get scale factors
  # Should we use the mean of the WT samples?
  #gmeans <- apply(norm_peaks[,6:(ncol(norm_peaks))], 1, function(x){exp(mean(log(x)))})
  gmeans <- apply(norm_peaks[,6:ncol(promoter_df)], 1, function(x){exp(mean(log(x)))})
  gnormed_counts <- norm_peaks[,6:(ncol(norm_peaks))] / gmeans
  scale_factors <- round(apply(gnormed_counts, 2, median),5)
  # adjust overlap_matrix by scale factors
  normalized_matrix <- round(t(t(overlap_matrix)/scale_factors),3)
  # Save the normalization peaks for later use
  write.table(norm_peaks, file = paste(dirname(output_filename), "/", "housekeeping_peaks_new.bed", sep=""), 
              sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
  # Also write the scale factors
  write.table(scale_factors, file = paste(dirname(output_filename), "/", "scale_factors_new.txt", sep=""),
              sep = '\t', col.names = FALSE, row.names = TRUE, quote = FALSE)
  
}else if(normalization == "peaks"){
  ### Reads-In-Peaks Normalization ###
  # Normalize all samples to have the same number of reads in peaks
  normalized_matrix <- cpm(overlap_matrix)
}

overlap_matrix <- cbind(as.character(peak_df$peak_name), overlap_matrix)
if(exists("normalized_matrix")){
  overlap_matrix <- cbind(as.character(peak_df$peak_name), normalized_matrix)
}

# Write counts matrix
write.table(overlap_matrix, file = output_filename, 
            sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
