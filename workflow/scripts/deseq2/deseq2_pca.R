#### PCA Plot to summarise differential expression analysis #### 

# Set working directory
setwd("~/OneDrive/Documents/Uni/III/Project/github/celegans-polyA/workflow")

# Load packages
library(DESeq2)
library(ggplot2)
library(stringr)

# Load rlog_dds object made previously
# load('output/deseq2_data/altadena/altadena_rlog_dds.RData')
load(snakemake@input$rlog_dds)

# Perform PCA analysis
pcaData <- DESeq2::plotPCA(rlog_dds, intgroup = "diet", returnData = TRUE)

# Calculate percentage variation explained by each of the principle components
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Legend

# Reorder
pcaData$diet <- factor(pcaData$diet, levels = c("as", "bp", "hb101", "op50", "pf", "m9"))

# Labels
legend_labels <- c(expression(italic("Acinetobacter schindleri")),
                   expression(italic("Bacillus pumilus")),
                   expression(paste(italic("Escherichia coli"), " strain HB101")),
                   expression(paste(italic("Escherichia coli"), " strain OP50")),
                   expression(italic("Pseudomonas fragi")),
                   "Starvation")

# Plot
# tiff('output/deseq2_plots/altadena/altadena_deseq2_pca.tiff', width = 1600, height = 800, res = 300)
tiff(snakemake@output$pca_plot, width = 1600, height = 800, res = 300)
ggplot2::ggplot(pcaData) +
  geom_point(mapping = aes(x = PC1, y = PC2, color = diet), size = 2) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle(stringr::str_to_title(snakemake@wildcards$location)) +
  labs(color = "Diet") +
  theme(text = element_text(size=10)) +
  scale_color_brewer(labels = legend_labels, palette = "Set2") +
  theme_light()
dev.off()
