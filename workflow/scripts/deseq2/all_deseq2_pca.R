#### PCA Plot to summarise differential expression analysis #### 

# Load packages
library(DESeq2)
library(ggplot2)

# Load rlog_dds object made previously
#vst_dds <- readRDS('output/deseq2_data/all/all_vst_dds.rds')
vst_dds <- readRDS(snakemake@input$vst_dds)

# Perform PCA analysis
pcaData <- DESeq2::plotPCA(vst_dds, intgroup = c('diet', 'location'), returnData = TRUE)

# Calculate percentage variation explained by each of the principle components
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Legend

# Reorder
pcaData$diet <- factor(pcaData$diet, levels = c("as", "bp", "hb101", "op50", "pf", "m9"))

pcaData$location <- factor(pcaData$location, levels = c('bristol', 'altadena'))

# Labels
diet_legend_labels <- c(expression(italic("Acinetobacter schindleri")),
                   expression(italic("Bacillus pumilus")),
                   expression(paste(italic("Escherichia coli"), " strain HB101")),
                   expression(paste(italic("Escherichia coli"), " strain OP50")),
                   expression(italic("Pseudomonas fragi")),
                   "Starvation")

location_legend_labels <- c('N2', 'PS2025')

# Plot
# tiff('output/deseq2_plots/all/all_deseq2_pca.tiff', width = 1600, height = 1000, res = 300)
ggplot2::ggplot(pcaData) +
  geom_point(mapping = aes(x = PC1, y = PC2, color = diet, shape = location), size = 2) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  guides(col=guide_legend("Diet"),
         shape=guide_legend("Strain")) +
  scale_color_brewer(labels = diet_legend_labels, palette = "Set2") +
  scale_shape_manual(labels = location_legend_labels, values = c(16, 18)) +
  theme_light() +
  theme(text = element_text(size = 10), panel.grid = element_blank()) +
  # guides(shape = FALSE) +
  ggtitle('All')

ggsave(filename = snakemake@output$pca_plot, device = 'tiff', width = 15,
       height = 7.5, units = 'cm', dpi = 320)
