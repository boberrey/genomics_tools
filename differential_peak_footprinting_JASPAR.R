# Differential peak footprinting

# (Ask some people if this is a crazy thing to do maybe...)

library(ggplot2)
library(ggrastr)
library(reshape2)
library(GenomicRanges)
library(GenomicAlignments)
library(Rsamtools)

source("/home/ben/git_clones/ATAC_tools/plotting_config.R")


footprint_file <- "output/footprinting/precomputed/per_sample_footprints.rds"
bias_file <- "output/footprinting/precomputed/bias_matrices.rds"
motif_kmer_file <- "/raid/USRdirs/ben/downloaded_genomic_data/mm10/jaspar_motif_hexamer_freq_matrices_3.rda"
diff_file <- "output/DESeq_diff_peaks/Hex15_diff_peaks.txt"
wt_bams <- Sys.glob("output/bams/deduped/WT_S_*.bam")
n_cores <- 5
g <- "mm10"

tf_flank <- 250

plot_dir <- "output/plots/TF_footprinting/JASPAR/Hex15_diff_peaks/"

# Get correct genome
if(g == "mm10"){
  library(BSgenome.Mmusculus.UCSC.mm10)
  gnome <- BSgenome.Mmusculus.UCSC.mm10
  sp <- "Mus musculus"
} else if(g == "hg38"){
  library(BSgenome.Hsapiens.UCSC.hg38)
  gnome <- BSgenome.Hsapiens.UCSC.hg38
  sp <- "Homo sapiens"
} else {
  print("Invalid genome. Quitting...")
  quit()
}

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

sp <- "Mus musculus"
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

# Specify parameters for reading in Bam files:
bam_param <- ScanBamParam(mapqFilter = 30)


# Function to get insertions
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
  reads <- GenomicRanges::shift(reads, ifelse(strand(reads) == '+', 4, -5))
  # Get overlaps
  return(reads)
}

insertion_ranges <- mclapply(wt_bams, function(x) get_insertions_from_bam(x, bam_param, gnome), mc.cores = n_cores)

wt_insertions <- do.call(c, insertion_ranges)


#-----------------------------------------------------
# Split motifs into those that fall into differential
# or constitutive peaks
#-----------------------------------------------------

# Now get motif positions in each peak set:
const_motif_matches <- matchMotifs(motifs, diff_const, genome = gnome, out = "positions")
diff_motif_matches <- matchMotifs(motifs, diff_down, genome = gnome, out = "positions")

# Resize motifs to indicated window:
resized_const_motifs <- endoapply(const_motif_matches, function(x) resize(x, width = (tf_flank*2 + 1), fix="center"))
resized_diff_motifs <- endoapply(diff_motif_matches, function(x) resize(x, width = (tf_flank*2 + 1), fix="center"))


split_by_strand <- function(mat, motif_ranges, motif_width, s){
  if(!(any(as.vector(strand(motif_ranges) == s)))){
    return(Matrix(0, nrow = 1, ncol = motif_width))
  }
  stranded <- as.matrix(mat[as.vector(strand(motif_ranges) == s),])
  if(dim(stranded)[2] == 1){
    return(t(stranded))
  }
  stranded
}

per_motif_footprint <- function(ins_cov, motif_gr, motif_width){
  # Get chromosomes to iterate over
  chroms <- seqlevels(motif_gr)
  motif_cov_list <- list()
  for(i in seq_along(chroms)){
    ch <- chroms[i]
    chr_motifs <- motif_gr[seqnames(motif_gr) == ch]
    # Skip chromosomes without any motifs
    if(length(chr_motifs) < 1){
      motif_cov_list[[i]] <- Matrix(0, nrow = motif_width, ncol = 1)
      next
    }
    m_view <- as.matrix(Views(ins_cov[ch], ranges(chr_motifs)))
    pos_view <- split_by_strand(m_view, chr_motifs, motif_width, '+')
    neg_view <- split_by_strand(m_view, chr_motifs, motif_width, '-')
    
    motif_cov_list[[i]] <- colSums(pos_view) + rev(colSums(neg_view))
  }
  Reduce(`+`, motif_cov_list)
}

# Calculate each motif footprint for each sample:
get_footprints <- function(ins_cov, motifs, n_cores=5){
  motif_width <- unique(sapply(motifs, function(x) unique(width(x))))
  # Iterate over motifs:
  motif_footprints <- mclapply(motifs, function(x) per_motif_footprint(ins_cov, x, motif_width), mc.cores = n_cores)
  #return(motif_footprints)
  footprints <- Reduce(cbind, motif_footprints)
  colnames(footprints) <- names(motifs)
  return(footprints)
}

# Calculate footprints in each group (slow)
wt_coverage <- coverage(wt_insertions)
const_footprints <- get_footprints(wt_coverage, resized_const_motifs, n_cores)
diff_footprints <- get_footprints(wt_coverage, resized_diff_motifs, n_cores)


#-----------------------------------------------------
# Proceed with plotting the two footprint sets
#-----------------------------------------------------


# Function to normalize insertion tracks to their flanking regions
normMat <- function(x, flank = 50){
  t(t(x) * (1/colMeans(x[c(1:flank, (nrow(x)-flank+1):nrow(x)),])))
}

# Function to get trimmed mean
trimmedMean <- function(x, trim_pct=0.1){
  lower_lim <- quantile(x, trim_pct)
  upper_lim <- quantile(x, 1-trim_pct)
  mean(x[(x > lower_lim) & (x < upper_lim)])
}

# Plotting function
# (Add something to show how many motifs in each group?)
plot_footprints <- function(footprint_mat, samp_labels, title="", ylabel=""){
  xpos <- seq(from=-floor(nrow(footprint_mat)/2), length.out = nrow(footprint_mat))
  footprint_mat <- as.data.frame(cbind(footprint_mat, xpos))
  # footprint and bias matrices shouls have sample names as column names
  melt_df <- melt(footprint_mat, id.vars = 'xpos')
  colnames(melt_df) <- c("xpos", "fp_name", "ins_signal")
  footp <- (
    ggplot(data = melt_df, aes(x=xpos, y=ins_signal, group=fp_name, color=fp_name))
    + geom_line()
    + theme_BOR()
    + xlab("Distance to Motif Center")
    + ylab(ylabel)
    + ggtitle(title)
    + scale_x_continuous(expand = c(0.005,0))
    #+ ylim(0.5, 1.5)
    + scale_colour_BOR(name = "Sample", labels=samp_labels)
    + theme(legend.justification = c(1, 1), 
            legend.position = c(1, 1), 
            legend.box.margin=margin(c(10,10,10,10)))
    + coord_cartesian(
      expand = FALSE, 
      ylim = c(0.8*quantile(melt_df$ins_signal,0.0001),1.2*quantile(melt_df$ins_signal,0.995)), 
      xlim = c(min(melt_df$xpos),max(melt_df$xpos))
    ) 
  )
  footp
}


# Read in necessary data
footprint_data <- list(const_footprints, diff_footprints)
names(footprint_data) <- c("const", "diff")
bias_data <- readRDS(bias_file)
# Warning: may be very large file ( ~16 Mb per matrix)
motif_kmer_data <- readRDS(motif_kmer_file)


bias_names <- c("WT_S", "WT_S")
samp_names <- c("const", "diff")


for(i in seq_along(names(motif_kmer_data))){
  motif_name <- names(motif_kmer_data)[i]
  kmer_freq <- motif_kmer_data[motif_name][[1]]
  
  # get the number of motifs in each group
  const_motifs <- length(resized_const_motifs[[motif_name]])
  diff_motifs <- length(resized_diff_motifs[[motif_name]])
  sample_labs <- c(paste(samp_names[1], as.character(const_motifs)), 
                   paste(samp_names[2], as.character(diff_motifs)))
  
  # Calculate expected insertions (kmer biased)
  bias_tracks <- Reduce("cbind", lapply(bias_data[bias_names], function(x) x[,3]))
  colnames(bias_tracks) <- samp_names
  bias_tracks <- t(kmer_freq) %*% bias_tracks
  exp_fp <- normMat(bias_tracks)
  
  # Retrieve actual insertions
  footprints <- Reduce("cbind", lapply(footprint_data[samp_names], function(x) x[,motif_name]))
  colnames(footprints) <- samp_names
  obs_fp <- normMat(footprints)
  
  # Observed / expected
  obsOverExp_fp <- obs_fp / exp_fp
  
  # Observed - expected
  obsMinusExp_fp <- obs_fp - exp_fp
  
  
  # Plot
  pdf(paste0(plot_dir, motif_name, "_diff.pdf"), width = 6, height = 5)
  print(plot_footprints(obs_fp, sample_labs, title = motif_name, ylabel = "Observed Insertions"))
  print(plot_footprints(exp_fp, sample_labs, title = motif_name, ylabel = "Expected Insertions"))
  print(plot_footprints(obsOverExp_fp, sample_labs, title = motif_name, ylabel = "Observed / Expected Insertions"))
  print(plot_footprints(obsMinusExp_fp, sample_labs, title = motif_name, ylabel = "Observed - Expected Insertions"))
  dev.off()
  # # Plot
  # colnames(exp_fp) <- sapply(samp_names, function(x) paste(x, "_bias"))
  # pdf(paste0(plot_dir, motif_name, "_narrow_new.pdf"), width = 6, height = 5)
  # print(plot_footprints(cbind(obs_fp[225:275,], exp_fp[225:275,]), title = motif_name, ylabel = "Observed Insertions"))
  # print(plot_footprints(exp_fp[225:275,], title = motif_name, ylabel = "Expected Insertions"))
  # print(plot_footprints(obsOverExp_fp[225:275,], title = motif_name, ylabel = "Observed / Expected Insertions"))
  # print(plot_footprints(obsMinusExp_fp[225:275,], title = motif_name, ylabel = "Observed - Expected Insertions"))
  # dev.off()
}



# Save Jeff's cross-footprint metrics
samp_list <- list()
for(i in seq_along(samp_names)){
  samp_list[[samp_names[i]]] <- list()
}

for(i in seq_along(names(motif_kmer_data))){
  motif_name <- names(motif_kmer_data)[i]
  kmer_freq <- motif_kmer_data[motif_name][[1]]
  
  # Calculate expected insertions (kmer biased)
  bias_tracks <- Reduce("cbind", lapply(bias_data[bias_names], function(x) x[,3]))
  colnames(bias_tracks) <- samp_names
  bias_tracks <- t(kmer_freq) %*% bias_tracks
  exp_fp <- normMat(bias_tracks)
  
  # Retrieve actual insertions
  footprints <- Reduce("cbind", lapply(footprint_data[samp_names], function(x) x[,motif_name]))
  colnames(footprints) <- samp_names
  obs_fp <- normMat(footprints)
  
  # Observed / expected
  obsOverExp_fp <- obs_fp / exp_fp
  
  # Observed - expected
  obsMinusExp_fp <- obs_fp - exp_fp
  
  # Calculate footprint metrics
  center_pos <- ceiling(nrow(exp_fp)/2)
  pwm_length <- ncol(as.matrix(motifs[motif_name][[1]]))
  base_flank <- floor(pwm_length/2) + 5
  base_positions <- (center_pos - base_flank):(center_pos + base_flank)
  
  flank_positions <- setdiff((center_pos-50):(center_pos+50), base_positions)
  
  backgroud_positions <- c(1:50, (nrow(exp_fp)-50):nrow(exp_fp))
  
  # Calculate for each sample:
  
  for(j in seq_along(samp_names)){
    samp <- samp_names[j]
    base_accessibility <- trimmedMean(obsOverExp_fp[base_positions,samp], trim_pct = 0.1)
    flank_accessibility <- mean(obsOverExp_fp[flank_positions,samp])
    bg_accessibility <- mean(obsOverExp_fp[backgroud_positions,samp])
    flank_fc <- log2(flank_accessibility/bg_accessibility)
    depth <- log2(base_accessibility/flank_accessibility)
    samp_list[[samp]][[motif_name]] <- c(bg_accessibility, flank_accessibility, base_accessibility, flank_fc, depth)
  }
}

# Convert into list of dataframes
result_dfs <- list()
for(i in seq_along(samp_names)){
  samp <- samp_names[i]
  rnames <- names(samp_list[[samp]])
  new_df <- as.data.frame(Reduce('rbind', samp_list[[samp]]))
  rownames(new_df) <- rnames
  colnames(new_df) <- c("bg_accessibility", "flank_accessibility", "base_accessibility", "flank_l2fc", "depth")
  result_dfs[[samp]] <- new_df
}

# Save footprint data
saveRDS(result_dfs, file = "output/footprinting/JASPAR/Hex15_diff_footprint_metrics.rds")