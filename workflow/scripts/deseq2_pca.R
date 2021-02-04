#### PCA Plot to summarise differential expression analysis #### 
library(DESeq2)
library(ggplot2)

# Load rlog_dds object made previously
#load('output/deseq2/bristol_rlog_dds.RData')
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
png(snakemake@output$pca_plot)
ggplot2::ggplot(pcaData) +
  geom_point(mapping = aes(x = PC1, y = PC2, color = diet), size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  labs(color = "Diet") +
  theme(text = element_text(size=14)) +
  scale_color_brewer(labels = legend_labels, palette = "Set2")
dev.off()