# run LOLA script

# LOLA analysis

# Fabian recommends using LOLA to find various associations with 
# available genomic datasets

# Usage:
# Rscript --vanilla run_LOLA_enrichment.R <diff_peak_file> <plot_dir> <data_dir> <species {'Hs', 'Mm'}>

suppressPackageStartupMessages({
  library(LOLA)
  library(GenomicRanges)
  library(ggplot2)
  library(reshape2)
})

# Get plotting settings
source("/home/ben/git_clones/ATAC_tools/plotting_config.R")

# Read in command line arguments
args <- commandArgs(trailingOnly = TRUE)
diff_file <- args[1]
plot_dir <- args[2]
data_dir <- args[3]
sp <- args[4] # Species {Mm, Hs}

# # args <- commandArgs(trailingOnly = TRUE)
# setwd("/raid/USRdirs/ben/phase_sep/hexanediol_ATAC/20190207_GM_and_HL_deep/")
# diff_file <- "output/diff_peaks/GM_15_diff_peaks.txt"
# plot_dir <- "output/plots/LOLA"
# data_dir <- "output/diff_peaks/LOLA_enrichments"
# sp <- "Hs" # Species {Mm, Hs}

if(sp == "Mm"){
  regionDB <- loadRegionDB("/raid/shr/Downloaded_data/regiondb/LOLACore/mm10/")
}else{
  regionDB <- loadRegionDB("/raid/shr/Downloaded_data/regiondb/LOLACore/hg38/")
}


diff_peaks_to_ranges <- function(diff_file){
  peak_df <- read.delim(diff_file)
  # Make GenomicRanges for peak df
  peak_ranges <- with(peak_df, GRanges(seqnames = chrom, 
                                       IRanges(start, end), 
                                       score = spm_score, 
                                       id = peak_name, 
                                       logFC = logFC, 
                                       logCPM = logCPM, 
                                       pvals = pvals, 
                                       fdr = fdr))
  return(peak_ranges)
}

diff_range <- diff_peaks_to_ranges(diff_file)

# Split range by groups (up, down)

fdr_cutoff <- 0.01
logFC_cutoff <- 0.5

get_diff_ranges <- function(diff_ranges, direction, fdr_cutoff=0.05, logFC_cutoff=0.5){
  result <- diff_ranges[(diff_ranges$fdr < fdr_cutoff) & 
                          (abs(diff_ranges$logFC) > logFC_cutoff)]
  if(direction == 'down'){
    return(result[result$logFC < 0,])
  }else{
    return(result[result$logFC > 0,])
  }
}

diff_down <- get_diff_ranges(diff_range, 'down', 
                              fdr_cutoff = fdr_cutoff, logFC_cutoff = logFC_cutoff)
diff_up <- get_diff_ranges(diff_range, 'up', 
                            fdr_cutoff = fdr_cutoff, logFC_cutoff = logFC_cutoff)


# Run LOLA
down_results <- runLOLA(diff_down, diff_range, regionDB = regionDB, cores=5)
up_results <- runLOLA(diff_up, diff_range, regionDB = regionDB, cores=5)


# Plot enrichments
plot_enrichments <- function(result, n=20){
  #valid_collections <- c("codex", "encode_tfbs", "cistrome_epigenome", "cistrome_cistrome", "ucsc_features")
  #plot_data <- result[result$collection %in% valid_collections][1:n]
  fixed_result <- result
  #max_p <- max(result$pValueLog[is.finite(result$pValueLog)])
  #fixed_result[which(!is.finite(fixed_result$pValueLog)),]$pValueLog <- max_p
  fixed_result[fixed_result$pValueLog > 300]$pValueLog <- 300
  plot_data <- fixed_result[1:n]
  plot_data$plot_name <- apply(plot_data, 1, function(x){
    s <- paste(paste(x[2], x[15], sep = "-"), paste(x[16], x[18], sep = "-"), sep='\n')
  })
  # Reorder factors
  plot_data$plot_name <- reorder(plot_data$plot_name, seq_along(plot_data$plot_name)[n:1])
  
  p <- (
    ggplot(data = plot_data, aes(x = plot_name, y = oddsRatio, fill = pValueLog))
    + geom_bar(stat="identity")
    + scale_fill_distiller(palette = "YlOrRd", direction = 1, limits=c(0, 300))
    + theme_BOR()
    + coord_flip()
    + ylab("Odds ratio")
    + xlab("")
    + theme(axis.text.y = element_text(hjust = 1, size = 7))
  )
  p
}

# Save plots
aspect_ratio <- 2
height <- 5

out_file_down <- paste(plot_dir, "/", sub('\\..*$', '', basename(diff_file)), "_down_LOLA.pdf", sep="")
down_p <- plot_enrichments(down_results)
ggsave(out_file_down, height = height, width = height*aspect_ratio)

out_file_up <- paste(plot_dir, "/", sub('\\..*$', '', basename(diff_file)), "_up_LOLA.pdf", sep="")
up_p <- plot_enrichments(up_results)
ggsave(out_file_up, height = height, width = height*aspect_ratio)

# Save data
out_file_down_data <- paste(data_dir, "/", sub('\\..*$', '', basename(diff_file)), "_down_LOLA.txt", sep="")
write.table(down_results, file = out_file_down_data, quote = F, sep='\t', row.names = F)
out_file_up_data <- paste(data_dir, "/", sub('\\..*$', '', basename(diff_file)), "_up_LOLA.txt", sep="")
write.table(up_results, file = out_file_up_data, quote = F, sep='\t', row.names = F)




