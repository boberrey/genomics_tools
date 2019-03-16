# Generate reads-in-peaks normalized bigwigs for ATAC data
# Based on the process described in Jeff's TCGA paper

# Usage: 
# generate_ATAC_signal_tracks.R <bam_file> <peaks_file> <output_filename> <genome> <bin_size> <normalization_style>
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomicAlignments)
  library(Rsamtools)
  library(rtracklayer)
  library(BSgenome)
})

# Read in command line arguments
args <- commandArgs(trailingOnly = TRUE)
bam_file <- BamFile(args[1])
peak_file <- args[2]
output_file <- args[3]
g <- args[4]
bin_size <- as.numeric(args[5])
norm_style <- args[6] # {'all', 'housekeeping', 'peaks'}
if(norm_style == 'housekeeping'){
  # If using 'housekeeping' normalization, provide the file containing those peaks
  housekeeping_file <- args[7]
}

# # troubleshooting
# bam_file <- BamFile("output/bams/deduped/mES_WT_1_S21.noMT.filtered.deduped.bam")
# peak_file <- "output/peaks/processed/filtered_fix_width_peaks.bed"
# output_file <- "output/coverage_data/tst.bw"
# g <- "mm10"
# bin_size <- as.numeric(100)
# norm_style <- 'housekeeping'
# if(norm_style == 'housekeeping'){
#   # If using 'housekeeping' normalization, provide the file containing those peaks
#   housekeeping_file <- "output/counts_matrix/scale_factors.txt"
# }

# Get correct genome
if(g == "mm10"){
  library(BSgenome.Mmusculus.UCSC.mm10)
  gnome <- BSgenome.Mmusculus.UCSC.mm10
} else if(g == "hg38"){
  library(BSgenome.Hsapiens.UCSC.hg38)
  gnome <- BSgenome.Hsapiens.UCSC.hg38
} else {
  print("Invalid genome. Quitting...")
  quit()
}



# Get tiled genome
chr_names <- seqnames(gnome)
chr_lengths <- seqlengths(gnome)
chr_ranges <- GRanges(seqnames = chr_names, 
                      IRanges(rep(1, length(chr_lengths)), chr_lengths),
                      seqlengths = chr_lengths)

# Remove non-standard chromosomes and Y and M
chr_ranges <- keepStandardChromosomes(chr_ranges, pruning.mode = 'tidy')
chr_ranges <- dropSeqlevels(chr_ranges, c("chrM", "chrY"), pruning.mode = 'tidy')

# Tile over genome (very large range produced)
tiled_chr_ranges <- tile(chr_ranges, width = bin_size)
names(tiled_chr_ranges) <- seqnames(chr_ranges)

# Get the merged and filtered peak set to get counts from
peak_df <- read.delim(peak_file, header=FALSE)
names(peak_df) <- c("chrom", "start", "end", "peak_name", "macs2_score", "spm_score")

# Make GenomicRanges for merged peaks
peak_ranges <- with(peak_df, GRanges(seqnames = chrom, 
                                     IRanges(start, end), 
                                     score = spm_score, 
                                     id = peak_name))


# Function to get shifted insertions
get_shifted_insertions_from_bam <- function(bam_file, bam_param=ScanBamParam(mapqFilter = 30)){
  # Read in bam file as GenomicRanges
  reads <- as(readGAlignments(bam_file, param = bam_param), "GRanges")
  # Resize to isolate the ends of fragments (the 5' end of each read)
  reads <- resize(reads, 1, fix = 'start')
  # Adjust insertion based on strand
  # ('+' stranded +4 bp, '-' stranded -5 bp)
  reads <- shift(reads, ifelse(strand(reads) == '+', 4, -5))
  # Remove chromosomes we dont want
  reads <- keepStandardChromosomes(reads, pruning.mode = 'tidy')
  reads <- dropSeqlevels(reads, c("chrM", "chrY"), pruning.mode = 'tidy')
  return(reads)
}

insertions <- get_shifted_insertions_from_bam(bam_file)


# Normalization
if(norm_style == 'housekeeping'){
  # load previously identified scale factors
  scale_factors <- read.delim(housekeeping_file, header=FALSE, row.names = 1)
  # Identify scale factor for this sample
  header <- basename(strsplit(bam_file$path, '\\.')[[1]][1])
  scale_factor <- scale_factors[header,1]
}else if(norm_style == 'peaks'){
  # Find scaling factor for reads-in-peaks normalization
  reads_in_peaks <- sum(countOverlaps(peak_ranges, insertions, ignore.strand=TRUE))
  # Previously they scaled to 30 million reads in peaks
  scale_factor <- reads_in_peaks/3e7
}else{
  # Or don't scale at all
  scale_factor <- 1
}



# Convert insertions range into 'coverage run-length encoding (Rle)'
ins_rle <- coverage(insertions)

# Get views of insertions over bined genome and calc sums
vSums <- viewSums(Views(ins_rle, as(tiled_chr_ranges, "RangesList")))

# Scale coverage by reads-in-peaks scaling factor
vSums <- vSums/scale_factor

# Re-order chromosomes to match tiled ranges
vSums <- vSums[names(tiled_chr_ranges)]

# Add vSums as score to tiled GRange
tiled_chr_ranges <- unlist(tiled_chr_ranges)
score(tiled_chr_ranges) <- unlist(vSums)
seqlengths(tiled_chr_ranges) <- chr_lengths[seqlevels(tiled_chr_ranges)]

# Save track in bigwig format
export.bw(tiled_chr_ranges, output_file)



