
library(ggplot2)
library(ggrastr)
library(reshape2)
source("/home/ben/git_clones/ATAC_tools/plotting_config.R")


footprint_file <- "output/footprinting/precomputed/per_sample_footprints_cisBP.rds"
bias_file <- "output/footprinting/precomputed/bias_matrices.rds"
motif_file <- "/raid/USRdirs/ben/downloaded_genomic_data/mm10/mouse_pwms_v2.rds"
motif_kmer_file <- "/raid/USRdirs/ben/downloaded_genomic_data/mm10/cisBP_motif_hexamer_freq_matrices.rds"
plot_dir <- "output/plots/TF_footprinting/cisBP/between_samples/"


# Get motifs:
motifs <- readRDS(motif_file)


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
plot_footprints <- function(footprint_mat, title="", ylabel=""){
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
    + scale_colour_BOR(name = "Sample")
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

# Function for plotting seq logo
plot_seq_logo <- function(TF_name, motifs, method="bits"){
  # First, get the motif ID
  mtx <- as.matrix(motifs[TF_name][[1]])
  # If any values < 0, assmume log-transformed PWM
  # (see https://github.com/GreenleafLab/chromVARmotifs)
  if(any(mtx < 0)){
    mtx <- exp(mtx)*0.25
  }
  # Fix ratio for different logo lengths
  ratio <- (1/ncol(mtx))*8
  p <- suppressMessages(
    ggplot()
    + geom_logo(mtx, method=method)
    + theme_logo()
    + ggtitle(TF_name)
    #+ coord_fixed(ratio = ratio)
  )
  p
}



# Read in necessary data
footprint_data <- readRDS(footprint_file)
bias_data <- readRDS(bias_file)
# Warning: may be very large file ( ~16 Mb per matrix)
motif_kmer_data <- readRDS(motif_kmer_file)


samp_names <- c("WT_S", "Hex60")


for(i in seq_along(names(motif_kmer_data))){
  motif_name <- names(motif_kmer_data)[i]
  kmer_freq <- motif_kmer_data[motif_name][[1]]
  
  # Calculate expected insertions (kmer biased)
  bias_tracks <- Reduce("cbind", lapply(bias_data[samp_names], function(x) x[,3]))
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
  
  # # Calculate footprint metrics
  # center_pos <- ceiling(nrow(exp_fp)/2)
  # pwm_length <- ncol(as.matrix(motifs[motif_name][[1]]))
  # base_flank <- floor(pwm_length/2) + 5
  # base_positions <- (center_pos - base_flank):(center_pos + base_flank)
  # 
  # flank_positions <- setdiff((center_pos-50):(center_pos+50), base_positions)
  # 
  # backgroud_positions <- c(1:50, (nrow(exp_fp)-50):nrow(exp_fp))
  # 
  # # Use trimmed mean for base, and regular mean for background and flank
  # base_accessibility <- trimmedMean(obsOverExp_fp[base_positions,])
  
  
  # Plot
  pdf(paste0(plot_dir, motif_name, ".pdf"), width = 6, height = 5)
  print(plot_footprints(obs_fp, title = motif_name, ylabel = "Observed Insertions"))
  print(plot_footprints(exp_fp, title = motif_name, ylabel = "Expected Insertions"))
  print(plot_footprints(obsOverExp_fp, title = motif_name, ylabel = "Observed / Expected Insertions"))
  print(plot_footprints(obsMinusExp_fp, title = motif_name, ylabel = "Observed - Expected Insertions"))
  sl1 <- plot_seq_logo(motif_name, motifs, method = "bits")
  sl2 <- plot_seq_logo(motif_name, motifs, method = "probability")
  grid.arrange(sl1, sl2, nrow=2)
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
samp_names <- names(bias_data)
samp_list <- list()
for(i in seq_along(samp_names)){
  samp_list[[samp_names[i]]] <- list()
}

for(i in seq_along(names(motif_kmer_data))){
  motif_name <- names(motif_kmer_data)[i]
  kmer_freq <- motif_kmer_data[motif_name][[1]]
  
  # Calculate expected insertions (kmer biased)
  bias_tracks <- Reduce("cbind", lapply(bias_data[samp_names], function(x) x[,3]))
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
saveRDS(result_dfs, file = "output/footprinting/cisBP/all_sample_footprint_metrics_cisBP.rds")
