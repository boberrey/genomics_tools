# Make fixed-width ATAC peaks and filter by score cutoff

# This script was implemented to perform Jeff Granja's method for obtaining 
# fixed-width ATAC peaks from a set of summits called on a single sample. Breifly:
# 
# 'For each sample, peak calling was performed on the Tn5-corrected single-base 
# insertions using the MACS2 callpeak command with parameters "--shift -75 --extsize 150
# --nomodel --call-summits --nolambda --keep-dup all -p 0.01". The peak summits 
# were then extended by 250 bp on either side to a final width of 501 bp, filtered
# by the ENCODE hg38 blacklist, and filtered to remove peaks that extend beyond
# the ends of chromosomes.
# 
# Overlapping peaks called within a single sample were handled using an iterative
# removal procedure. First, the most significant peak is kept and any peak that 
# directly overlaps with that significant peak is removed. Then, this process iterates
# to the next most significant peak and so on until all peaks have either been kept
# or removed due to direct overlap with a more significant peak...'


# Usage:
# get_fixed_width_peaks.R <peak_path> <out_path> <genome> <spm_cutoff> <extend_size> <n_cores>

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomicAlignments)
  library(preprocessCore)
})

# Read in command line arguments
args <- commandArgs(trailingOnly = TRUE)
peak_path <- args[1]
out_path <- args[2]
gnome_str <- args[3]
spm_cutoff <- as.numeric(args[4])
ext_size <- as.numeric(args[5])
n_cores <- as.numeric(args[6])

peak_files <- list.files(path = peak_path, pattern = "*summits.bed$", full.names = TRUE)


# Read in peak files and expand window
get_peak_df <- function(peak_file, window_size){
  peak_df <- read.delim(peak_file, header=FALSE)
  names(peak_df) <- c("chrom", "start", "end", "peak_name", "macs2_score")
  # Expand by window size
  peak_df$start <- peak_df$start - window_size
  peak_df$end <- peak_df$end + window_size
  peak_df
}

peak_dfs <- mclapply(peak_files, function(x) get_peak_df(x, ext_size), mc.cores=n_cores)


# Get genome information
gnome <- Seqinfo(genome = gnome_str)


# Remove overlapping peaks for each sample 
per_sample_peak_merge <- function(peak_df, spm_cutoff, gnome){
  # first convert peak set into GenomicRange
  peak_df <- peak_df[order(peak_df$macs2_score, decreasing = TRUE),]
  peak_ranges <- with(peak_df, GRanges(seqnames = chrom, 
                                       IRanges(start, end), 
                                       macs2_score = macs2_score,
                                       id = peak_name))
  # Remove peaks that extend past chromosome ends
  suppressWarnings(seqinfo(peak_ranges) <- gnome[seqlevels(peak_ranges)])
  idx <- GenomicRanges:::get_out_of_bound_index(peak_ranges)
  if(length(idx) > 0){
    peak_ranges <- peak_ranges[-idx]
  }

  # Iteratively remove peaks that overlap with the highest scoring peak
  reduced_ranges <- reduce(peak_ranges, ignore.strand=TRUE, min.gapwidth=0)
  while(length(reduced_ranges) != length(peak_ranges)){
    segments <- reduced_ranges[!reduced_ranges %in% peak_ranges]
    highest <- peak_ranges[findOverlaps(segments, peak_ranges, select = "first")]
    peak_ranges <- append(highest, peak_ranges[!peak_ranges %over% highest])
    reduced_ranges <- reduce(peak_ranges, ignore.strand=TRUE, min.gapwidth=0)
  }
  
  # Convert back to data frame
  peak_ranges <- sort(peak_ranges)
  new_df <- DataFrame(
    chrom=seqnames(peak_ranges),
    start=start(peak_ranges),
    end=end(peak_ranges),
    peak_name=peak_ranges$id,
    macs2_score=peak_ranges$macs2_score
  )
  # Calculate 'score per million' and filter by cutoff
  new_df$spm <- round(new_df$macs2_score / (sum(new_df$macs2_score)/1e6),5)
  new_df[new_df$spm > spm_cutoff,]
}


# Merge peaks for each sample and filter by spm:
filtered_dfs <- mclapply(peak_dfs, function(x) per_sample_peak_merge(x, spm_cutoff, gnome), mc.cores = n_cores)


# Save individual filtered peak files
generate_outfiles <- function(x) paste(out_path, "/", substr(basename(x), 1, regexpr("summits.bed", basename(x))[1]-1), "fix_width.bed", sep = "")
new_filenames <- unlist(lapply(peak_files, generate_outfiles))

for(i in 1:length(new_filenames)){
  filename <- new_filenames[i]
  df <- filtered_dfs[i]
  write.table(df, filename, quote=FALSE, sep='\t', row.names = FALSE, col.names = FALSE)
}


# Create merged peakset
all_sample_peak_merge <- function(peak_dfs){
  # Combine dataframes
  merged_df <- do.call("rbind", peak_dfs)
  # sort by spm and convert peak set into GenomicRange
  merged_df <- merged_df[order(merged_df$spm, decreasing = TRUE),]
  peak_ranges <- with(merged_df, GRanges(seqnames = chrom, 
                                       IRanges(start, end), 
                                       macs2_score = macs2_score,
                                       spm = spm,
                                       id = peak_name))
  
  # Get only unique peaks
  peak_ranges <- unique(peak_ranges)
  
  # Iteratively remove peaks that overlap with the highest scoring peak
  reduced_ranges <- reduce(peak_ranges, ignore.strand=TRUE, min.gapwidth=0)
  while(length(reduced_ranges) != length(peak_ranges)){
    segments <- reduced_ranges[!reduced_ranges %in% peak_ranges]
    highest <- peak_ranges[findOverlaps(segments, peak_ranges, select = "first")]
    peak_ranges <- append(highest, peak_ranges[!peak_ranges %over% highest])
    reduced_ranges <- reduce(peak_ranges, ignore.strand=TRUE, min.gapwidth=0)
  }
  # Convert back to data frame
  peak_ranges <- sort(peak_ranges)
  DataFrame(
    chrom=seqnames(peak_ranges),
    start=start(peak_ranges),
    end=end(peak_ranges),
    peak_name=peak_ranges$id,
    macs2_score=peak_ranges$macs2_score,
    spm=peak_ranges$spm
  )
}


merged_df <- all_sample_peak_merge(filtered_dfs)

# Save the merged peakset
write.table(merged_df, paste(out_path, "/merged_fix_width_peaks.bed", sep = ""), quote=FALSE, sep='\t', row.names = FALSE, col.names = FALSE)

