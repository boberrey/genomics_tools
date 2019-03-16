# Generate footprints for all samples accross all motifs

# Usage: 
# TF_footprint_all_samples.R <peak_file> <bam_dir> <output_dir> <tf_flank> <genome> <n_cores>

# In most cases use tf_flank = 250

suppressPackageStartupMessages({
  library(BSgenome)
  library(Biostrings)
  library(GenomicRanges)
  library(GenomicAlignments)
  library(Rsamtools)
  library(motifmatchr)
  library(JASPAR2018)
})

# Read in command line arguments (DEBUG)
args <- commandArgs(trailingOnly = TRUE)
peak_file <- args[1]
bam_dir <- args[2]
out_dir <- args[3]
tf_flank <- as.numeric(args[4]) # Usually 250
g <- args[5] # {'mm10', 'hg38'}
n_cores <- as.numeric(args[6])

# # Read in command line arguments (DEBUG)
# setwd("/raid/USRdirs/ben/phase_sep/hexanediol_ATAC/20190213_mESC_retry/")
# peak_file <- "output/peaks/processed/filtered_fix_width_peaks.bed"
# bam_dir <- "output/bams/deduped/"
# out_dir <- "output/footprinting/precomputed/"
# tf_flank <- 250
# g <- "mm10"
# n_cores <- as.numeric('10')

# Make directory if it doesn't already exist
dir.create(file.path(getwd(), out_dir), showWarnings = FALSE)

# Get correct genome
if(g == "mm10"){
  library(org.Mm.eg.db)
  library(BSgenome.Mmusculus.UCSC.mm10)
  gnome <- BSgenome.Mmusculus.UCSC.mm10
  sp <- "Mus musculus"
} else if(g == "hg38"){
  library(org.Hs.eg.db)
  library(BSgenome.Hsapiens.UCSC.hg38)
  gnome <- BSgenome.Hsapiens.UCSC.hg38
  sp <- "Homo sapiens"
} else {
  print("Invalid genome. Quitting...")
  quit()
}

print("Calculating genome-wide kmer frequency...")
# Create BSParams object containing the bsapply conditions for counting hexamer frequencies on standard chromosomes:
params <- new("BSParams", X = gnome, FUN = function(x){oligonucleotideFrequency(x, width=6)}, exclude = c("_", "M", "Y"))
# Run count and sum chromosome results
counts <- rowSums(as.data.frame(bsapply(params)))
genome_freqs <- counts / sum(counts)

print("Retrieving JASPAR motifs...")
# Get motifs:
getJasparMotifs <- function(species = "Homo sapiens",
                            collection = "CORE", ...) {
  opts <- list()
  opts["species"] <- species
  opts["collection"] <- collection
  opts <- c(opts, list(...))
  out <- TFBSTools::getMatrixSet(JASPAR2018::JASPAR2018, opts)
  if (!isTRUE(all.equal(TFBSTools::name(out), names(out))))
    names(out) <- paste(names(out), TFBSTools::name(out), sep = "_")
  return(out)
}

motifs <- getJasparMotifs(species = sp)


# Get the merged and filtered peak set
peak_df <- read.delim(peak_file, header=FALSE)

names(peak_df) <- c("chrom", "start", "end", "peak_name", "macs2_score", "spm_score")

# Make GenomicRanges for merged peaks
peak_ranges <- with(peak_df, GRanges(seqnames = chrom, 
                                     IRanges(start, end), 
                                     score = spm_score, 
                                     id = peak_name))

# Locate bam files
bam_files <- list.files(path = bam_dir, pattern = "\\.bam$", full.names = TRUE)

# Prepare for use with bioconductor I/O tools
headers <- sapply(strsplit(bam_files, '\\.'), function(l) l[[1]])
headers <- basename(headers)
bam_file_list <- lapply(bam_files, BamFile)
names(bam_file_list) <- headers

# Specify parameters for reading in Bam files:
bam_param <- ScanBamParam(mapqFilter = 30)


# Function to get counts in peaks
get_insertions_from_bam <- function(bam_file, bam_param, gnome){
  # Read in bam file as GenomicRanges
  reads <- as(readGAlignments(bam_file, param = bam_param), "GRanges")
  seqinfo(reads) <- seqinfo(gnome)[seqlevels(reads)]
  # Remove non-standard chromosomes and Y and M
  reads <- keepStandardChromosomes(reads, pruning.mode = 'tidy')
  reads <- dropSeqlevels(reads, c("chrM", "chrY"), pruning.mode = 'tidy')
  # Resize to isolate the ends of fragments (the 5' end of each read)
  reads <- resize(reads, 1, fix = 'start')
  # Adjust insertion based on strand
  # ('+' stranded +4 bp, '-' stranded -5 bp)
  reads <- shift(reads, ifelse(strand(reads) == '+', 4, -5))
  # Get overlaps
  return(reads)
}

print("Reading and merging insertions...")
insertion_ranges <- mclapply(bam_file_list, function(x) get_insertions_from_bam(x, bam_param, gnome), mc.cores = n_cores)

# Get individual group names
group_names <- sapply(names(insertion_ranges), function(x){
  # Assumes that technical replicates have {sample_identifier}_{tech rep #}_{seq samp #}
  # Split off the last two and get unique sample identifiers
  s <- strsplit(x, "_")
  paste(s[[1]][1:(length(s[[1]])-2)], collapse = "_")
})
group_names <- unique(group_names)

# Merge technical replicates before calculating bias and footprints:
merged_insertions <- list()
for(i in seq_along(group_names)){
  grp <- group_names[i]
  ingroup <- insertion_ranges[grep(paste('^', grp, sep=""), names(insertion_ranges))]
  merged_insertions[[grp]] <- Reduce("c", ingroup)
}


# Generate sample-specific Tn5 insertion bias matrix
get_insertion_bias <- function(ins_range, gnome, kmer_exp){
  ##################################################################
  # 'ins_range' is the Tn5-adjusted insertion sites from Bam file
  # 'gnome' is the appropriate BSgenome object
  # 'kmer_exp' is the gnome-calculated n-mer frequencies 
  ##################################################################
  # Infer kmer length from expected frequencies:
  k <- unique(sapply(names(kmer_exp), nchar))
  # Adjust to kmer window
  # (Jeff doesn't shift insertions based on strand, so I guess let's not either?)
  #reads <- trim(shift(resize(ins_range, width = k, fix = "center"), ifelse(strand(ins_range) == "+", 0, 1)))
  reads <- trim((resize(ins_range, width = k, fix = "center")))
  # Count kmer frequencies
  kmer_obs <- oligonucleotideFrequency(getSeq(gnome, reads), width=k, simplify.as="collapse")
  kmer_obs <- kmer_obs / sum(kmer_obs)
  kmer_OE <- kmer_obs[names(kmer_exp)] / kmer_exp
  return(cbind(kmer_obs, kmer_exp, kmer_OE))
}

print("Calculating bias matrices for all samples...")
# Calculate bias matrices for WT and experimental value (slow)
bias_matrices <- mclapply(merged_insertions, function(x) get_insertion_bias(x, gnome, genome_freqs), mc.cores = n_cores)

# Now get motif positions for all samples:
motif_matches <- matchMotifs(motifs, peak_ranges, genome = gnome, out = "positions")

# Resize motifs to indicated window:
resized_motifs <- endoapply(motif_matches, function(x) resize(x, width = (tf_flank*2 + 1), fix="center"))

# Calculate each motif footprint for each sample:
get_footprints <- function(insertions, motifs){
  # Get coverage:
  ins_cov <- coverage(insertions)
  # Iterate over motifs:
  motif_footprints <- list()
  for(i in seq_along(motifs)){
    # Get motif:
    motif_gr <- motifs[[i]]
    # Get chromosomes to iterate over
    chroms <- seqlevels(motif_gr)
    motif_cov_list <- list()
    for(j in seq_along(chroms)){
      ch <- chroms[j]
      chr_motifs <- motif_gr[seqnames(motif_gr) == ch]
      m_view <- as.matrix(Views(ins_cov[ch], ranges(chr_motifs)))
      pos_view <- m_view[as.vector(strand(chr_motifs) == "+"),]
      neg_view <- m_view[as.vector(strand(chr_motifs) == "-"),]
      
      motif_cov_list[[j]] <- colSums(pos_view) + colSums(neg_view[,ncol(neg_view):1])
    }
    motif_footprints[[i]] <- Reduce(`+`, motif_cov_list)
  }
  footprints <- Reduce(cbind, motif_footprints)
  colnames(footprints) <- names(motifs)
  return(footprints)
}

print("Calculating all footprints...")
footprints <- mclapply(merged_insertions, function(x) get_footprints(x, resized_motifs), mc.cores = n_cores)

# Save all outputs:
bias_file <- paste(out_dir, "/", "bias_matrices.rds",sep="")
saveRDS(bias_matrices, file=bias_file)

footprint_file <- paste(out_dir, "/", "per_sample_footprints.rds",sep="")
saveRDS(footprints, file=footprint_file)