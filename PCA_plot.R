
library(ggplot2)


# Do a quick PCA plot?
count_matrix <- read.delim("output/counts_matrix/raw_counts_table.txt")
rownames(count_matrix) <- count_matrix[,1]
count_matrix <- count_matrix[,c(-1, -6, -7, -8, -9)]

input_matrix <- normalize.quantiles(cpm(count_matrix, log = TRUE, prior.count = 3))

colnames(input_matrix) <- colnames(count_matrix)
rownames(input_matrix) <- rownames(count_matrix)

pc <- prcomp(t(input_matrix))
variance <- pc$sdev^2 / sum(pc$sdev^2)

var_data <- data.frame(PC = c(1:length(variance)), variance)

var_plot <- (
  ggplot(var_data, aes(x = PC, y = variance))
  + geom_bar(stat = 'identity')
  + theme_bw()
  + theme(text = element_text(size=12),
          panel.grid.major=element_blank(), 
          panel.grid.minor= element_blank(), 
          aspect.ratio = 1)
  + ylab("Fraction of variance")
  + scale_x_discrete(name = "PC", limits = var_data$PC)
)
var_plot
ggsave("pca_var_explained_subset.pdf")

pc_data <- cbind(data.frame(Sample=rownames(pc$x)), pc$x)
pc_data$group <- c("1.5% Hex-1,6", "1.5% Hex-1,6", 
                   "1.5% Hex-2,5", "1.5% Hex-2,5", 
                   #"5% Hex-1,6", "5% Hex-1,6", 
                   #"5% Hex-2,5", "5% Hex-2,5", 
                   "WT", "WT")

max_pc <- 4

color_column <- NA
shape_column <- NA

for(i in seq(1, max_pc - 1)) {
  for(j in seq(i + 1, max_pc)) {
    pcplot = (ggplot(pc_data, aes_string(x = paste0("PC",i), y = paste0("PC",j), group = "group", color="group")) 
              + geom_point() 
              + theme_bw() 
              + theme(text = element_text(size=12),
                      panel.grid.major=element_blank(), 
                      panel.grid.minor= element_blank()) 
              + scale_color_brewer(palette="Set1") 
              + coord_equal()
              )
    #if(color_column == "None") pcplot = pcplot + scale_color_manual(values="#756bb1", guide = FALSE)
    #if(shape_column == "None") pcplot = pcplot + scale_shape_discrete(guide = FALSE)
    ggsave(paste0("pc", i, "vs", j, "_subset.pdf"))
  }
}

