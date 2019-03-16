
library(JASPAR2016)
library(motifmatchr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(SummarizedExperiment)
library(Matrix)



# Get JASPAR motifs using function from Alica's ChromVAR
getJasparMotifs <- function(species = "Homo sapiens", 
                            collection = "CORE", ...) {
  opts <- list()
  opts["species"] <- species
  opts["collection"] <- collection
  opts <- c(opts, list(...))
  out <- TFBSTools::getMatrixSet(JASPAR2016::JASPAR2016, opts)
  if (!isTRUE(all.equal(TFBSTools::name(out), names(out)))) 
    names(out) <- paste(names(out), TFBSTools::name(out), sep = "_")
  return(out)
}

motifs <- getJasparMotifs(species = "Mus musculus")


diff_peaks <- "hexanediol_16_5pct_diff_peaks.txt"
peak_df <- read.delim(diff_peaks)

# Make GenomicRanges for peak df
peak_ranges <- with(peak_df, GRanges(seqnames = chrom, 
                                     IRanges(start, end), 
                                     score = spm_score, 
                                     id = peak_name, 
                                     logFC = logFC, 
                                     logCPM = logCPM, 
                                     pvals = pvals, 
                                     fdr = fdr))


# Split range by groups (constitutive, up, down)

fdr_cutoff <- 0.05
logFC_cutoff <- 1

diff_ranges <- peak_ranges[(peak_ranges$fdr < fdr_cutoff) & 
                             (abs(peak_ranges$logFC) > logFC_cutoff)]

diff_up_ranges <- diff_ranges[diff_ranges$logFC > 0,]
diff_down_ranges <- diff_ranges[diff_ranges$logFC < 0,]


const_ranges <- peak_ranges[(abs(peak_ranges$logFC) < 1),]


# Function to get motif counts across whole range
total_motif_counts <- function(ranges, motifs, genome=BSgenome.Mmusculus.UCSC.mm10){
  # identify motif matches
  motif_ix <- matchMotifs(motifs, ranges, genome = genome)
  sums <- apply(assay(motif_ix), 2, sum)
  names(sums) <- motif_ix$name
  sums
}

# total_counts <- total_motif_counts(peak_ranges, motifs)
# const_counts <- total_motif_counts(const_ranges, motifs)
# up_counts <- total_motif_counts(diff_up_ranges, motifs)
# down_counts <- total_motif_counts(diff_down_ranges, motifs)

motif_df <- DataFrame(down_counts, total_counts)

# Perform hypergeometric tests for each motif
motif_hypergeom_tests <- function(range_1, range_2, motifs){
  n1 <- length(range_1)
  n2 <- length(range_2)
  counts_1 <- total_motif_counts(range_1, motifs)
  counts_2 <- total_motif_counts(range_2, motifs)
  df <- DataFrame(counts_1, counts_2)
  pvals <- round(apply(df, 1, function(x) motif_hypergeom(x, n1, n2)), 3)
  log2_FC <- round(log2((counts_1/n1)/(counts_2/n2)), 3)
  cbind(df, pvals, log2_FC)
}

motif_hypergeom <- function(row, n1, n2){
  m1 <- row[1] # Number of motifs in set
  m2 <- row[2] # Total number of motifs
  -phyper(m1 - 1, m2, n2, n1, log.p = T, lower.tail = FALSE)
}

diff_motifs_down <- motif_hypergeom_tests(diff_down_ranges, peak_ranges, motifs)
diff_motifs_up <- motif_hypergeom_tests(diff_up_ranges, peak_ranges, motifs)

tail(tst[order(tst$pvals, decreasing = T),], 10)


