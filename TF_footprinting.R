# More motif analysis (and footprinting?)

library(BSgenome)
library(Biostrings)
library(GenomicRanges)
library(GenomicAlignments)
library(Rsamtools)
library(motifmatchr)
library(JASPAR2018)
source("/home/ben/git_clones/ATAC_tools/plotting_config.R")

# # Read in command line arguments
# args <- commandArgs(trailingOnly = TRUE)
# peak_file <- args[1]
# bam_file_path <- args[2]
# output_filename <- args[3]
# n_cores <- as.numeric(args[4])
# g <- args[5] # {'mm10', 'hg38'}

# Read in command line arguments (DEBUG)
setwd("/raid/USRdirs/ben/phase_sep/hexanediol_ATAC/20190213_mESC_retry/")
diff_file <- "output/DESeq_diff_peaks/Hex15_diff_peaks.txt"
wt_bams <- Sys.glob("output/bams/deduped/WT_S_*.bam")
exp_bams <- Sys.glob("output/bams/deduped/Hex15*.bam")
motif_bias_file <- "/raid/USRdirs/ben/downloaded_genomic_data/mm10/jaspar_motif_hexamer_freq_matrices.rda"
plot_dir <- "output/plots/TF_footprinting"
n_cores <- as.numeric('4')
tf_flank <- 250
g <- "mm10"

# Make directory if it doesn't already exist
dir.create(file.path(getwd(), plot_dir), showWarnings = FALSE)

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

# Create BSParams object containing the bsapply conditions for counting hexamer frequencies on standard chromosomes:
params <- new("BSParams", X = gnome, FUN = function(x){oligonucleotideFrequency(x, width=6)}, exclude = c("_", "M", "Y"))
# Run count and sum chromosome results
counts <- rowSums(as.data.frame(bsapply(params)))
genome_freqs <- counts / sum(counts)


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


# Get the merged and filtered peak set to get counts from
diff_df <- read.delim(diff_file, header=T)

# Differential peaks first:
diff_peaks_to_ranges <- function(diff_file){
  peak_df <- read.delim(diff_file)
  # Make GenomicRanges for peak df
  peak_ranges <- with(peak_df, GRanges(seqnames = chrom, 
                                       IRanges(start, end), 
                                       gene_symbol = as.character(symbol),
                                       tss_dist = distance,
                                       id = as.character(peak_name), 
                                       logFC = log2FC, 
                                       logCPM = log2CPM, 
                                       fdr = padj))
  return(peak_ranges)
}
diff_range <- diff_peaks_to_ranges(diff_file)

fdr_cutoff <- 0.01
fold_change <- -0.5
diff_down <- diff_range[(diff_range$fdr < fdr_cutoff) & (diff_range$logFC < fold_change)]
diff_const <- setdiff(diff_range, diff_down)

# Prepare for use with bioconductor I/O tools

wt_bams <- lapply(wt_bams, BamFile)
exp_bams <- lapply(exp_bams, BamFile)

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

insertion_ranges <- mclapply(c(wt_bams, exp_bams), function(x) get_insertions_from_bam(x, bam_param, gnome), mc.cores = n_cores)

wt_insertions <- do.call(c, insertion_ranges[1:length(wt_bams)])
exp_insertions <- do.call(c, insertion_ranges[(1+length(wt_bams)):length(insertion_ranges)])

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
  reads <- trim(shift(resize(ins_range, width = k, fix = "center"), ifelse(strand(ins_range) == "+", 0, 1)))
  # Count kmer frequencies
  kmer_obs <- oligonucleotideFrequency(getSeq(gnome, reads), width=k, simplify.as="collapse")
  kmer_obs <- kmer_obs / sum(kmer_obs)
  kmer_OE <- kmer_obs / kmer_exp
  return(cbind(kmer_obs, kmer_exp, kmer_OE))
}

# Calculate bias matrices for WT and experimental value (slow)
wt_bias_mat <- get_insertion_bias(wt_insertions, gnome, genome_freqs)
exp_bias_mat <- get_insertion_bias(exp_insertions, gnome, genome_freqs)

# Read in pre-computed motif kmer bias tables
motif_kmer_bias <- readRDS(motif_bias_file)

# Get motif positions
motif_matches <- matchMotifs(motifs, diff_range, genome = gnome, out = "positions")

resized_motifs <- endoapply(motif_matches, function(x) resize(x, width = (tf_flank*2 + 1), fix="center"))

# Calculate motif coverage


wt_cov <- coverage(wt_insertions)
exp_cov <- coverage(exp_insertions)

# Split wt insertions into differential vs not?
wt_const_ins <- subsetByOverlaps(wt_insertions, diff_const, maxgap = 100)
wt_diff_ins <- subsetByOverlaps(wt_insertions, diff_down, maxgap = 100)
wt_const_ins <- sort(sample(wt_const_ins, size = length(wt_diff_ins), replace = F))
wt_const_cov <- coverage(wt_const_ins)
wt_diff_cov <- coverage(wt_diff_ins)

normMat <- function(x, flank = 50){
  # Footprint matrix normalization from Jeff
  t(t(x) * (1/colMeans(x[c(1:flank, (nrow(x)-flank+1):nrow(x)),])))
}

get_footprint <- function(ins_cov, motif_gr){
  # Get chromosomes to iterate over
  chroms <- seqlevels(motif_gr)
  motif_cov_list <- list()
  for(i in seq_along(chroms)){
    ch <- chroms[i]
    chr_motifs <- motif_gr[seqnames(motif_gr) == ch]
    #m_view <- as.matrix(BSgenome::Views(ins_cov[ch], start = start(chr_motifs), end = end(chr_motifs)))
    m_view <- as.matrix(Views(ins_cov[ch], ranges(chr_motifs)))
    pos_view <- m_view[as.vector(strand(chr_motifs) == "+"),]
    neg_view <- m_view[as.vector(strand(chr_motifs) == "-"),]
    
    motif_cov_list[[i]] <- colSums(pos_view) + colSums(neg_view[,ncol(neg_view):1])
  }
  return(Reduce(`+`, motif_cov_list))
}

motif_name <- "MA0514.1_Sox3"
tst_motif <- unlist(resized_motifs[motif_name])
# motif_cov_wt <- get_footprint(wt_cov, tst_motif)
# motif_cov_exp <- get_footprint(exp_cov, tst_motif)

motif_cov_wt <- get_footprint(wt_const_cov, tst_motif)
motif_cov_exp <- get_footprint(wt_diff_cov, tst_motif)

mt <- normMat(cbind(motif_cov_wt, motif_cov_exp))

p <- (
  ggplot()
  + geom_line(aes(x=(-250:250), y=mt[,1]), color="red")
  + geom_line(aes(x=(-250:250), y=mt[,2]), color="blue")
  + theme_BOR()
  + xlab("Distance to Motif Center")
  + ylab("Normalized Insertions")
  + scale_x_continuous(expand = c(0.005,0))
)
p

