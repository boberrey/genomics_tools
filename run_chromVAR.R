# Run chromVAR

library(chromVAR)
library(SummarizedExperiment)
library(motifmatchr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BiocParallel)
library(JASPAR2016)

# Prepare parallelization
register(MulticoreParam(10))


peak_file <- "output/peaks/processed/filtered_fix_width_peaks.bed"
peaks <- getPeaks(peak_file, sort_peaks = TRUE)

# How to use chromVAR for bulk data?

bam_file_path <- "output/bams/deduped"
bam_files <- unlist(list.files(path = bam_file_path, pattern = "\\.bam$", full.names = TRUE))

fragment_counts <- getCounts(bam_files, peaks, 
                             paired = TRUE, by_rg = TRUE, 
                             colData = DataFrame(treatment = c("1.5% Hex-1,6", "1.5% Hex-1,6", 
                                                               "1.5% Hex-2,5", "1.5% Hex-2,5", 
                                                               "5% Hex-1,6", "5% Hex-1,6", 
                                                               "5% Hex-2,5", "5% Hex-2,5", 
                                                               "WT", "WT")))


fragment_counts <- addGCBias(fragment_counts, genome = BSgenome.Mmusculus.UCSC.mm10)

# The fraction of reads in peaks is abysmal for the high treatment samples. 
filtering_plot <- filterSamplesPlot(fragment_counts, min_depth = 1500, min_in_peaks = 0.05, use_plotly = FALSE)
filtering_plot

motifs <- getJasparMotifs(species = "Mus musculus")

motif_ix <- matchMotifs(motifs, fragment_counts, genome = BSgenome.Mmusculus.UCSC.mm10)


# Compute deviations
dev <- computeDeviations(object = fragment_counts, annotations = motif_ix)


variability <- computeVariability(dev)
plotVariability(variability, use_plotly = FALSE)


# tSNE

tsne_results <- deviationsTsne(dev, threshold = 1.5, perplexity = 3)
tsne_plots <- plotDeviationsTsne(dev, tsne_results, 
                                 sample_column = "treatment", 
                                 annotation_name =c("Pou5f1::Sox2", "Sox2", "Bach1::Mafk", "Nobox", "Neurog1", "Lhx8", "RUNX1", "Dlx1"), 
                                 shiny = FALSE)



