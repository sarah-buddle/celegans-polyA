#### Plot percentage expressed genes across all samples

# Set working directory
setwd("~/OneDrive/Documents/Uni/III/Project/github/celegans-polyA/workflow")

# Load packages
library(ggplot2)

# Load data objects
#load('output/expressed_genes/altadena/altadena_expressed_genes.RData')
load(snakemake@input$full_dds_altadena)
expressed_genes1 <- expressed_genes

#load('output/expressed_genes/bristol/bristol_expressed_genes.RData')
load(snakemake@input$full_dds_bristol)
expressed_genes2 <- expressed_genes

# Merge into single data frame

expressed_genes_all <- rbind(expressed_genes1, expressed_genes2)
expressed_genes_all <- aggregate(expressed_genes_all,
                                 by = list(expressed_genes_all$location,
                                          expressed_genes_all$diet),
                                 FUN = mean)
expressed_genes_all <- expressed_genes_all[, c('Group.1', 'Group.2', 
                                               'percentage_expressed_genes')]
colnames(expressed_genes_all) = c('location', 'diet', 'percentage_expressed_genes')

# Plot percentages

# Legend
# Reorder
expressed_genes_all$diet <- factor(expressed_genes_all$diet,
                               levels = rev(c("as", "bp", "hb101", "op50", "pf", "m9")))
# labels
x_labels <- rev(c(expression(italic("Acinetobacter schindleri")),
                   expression(italic("Bacillus pumilus")),
                   expression(paste(italic("Escherichia coli"), " strain HB101")),
                   expression(paste(italic("Escherichia coli"), " strain OP50")),
                   expression(italic("Pseudomonas fragi")),
                   "Starvation"))

legend_labels <- c('Altadena', 'Bristol')


# Plot
#tiff('output/deseq2_plots/all/percentage_plot_all.tiff', width = 1600, height = 800, res = 300)
tiff(snakemake@output$percentage_plot, width = 1600, height = 800, res = 300)
ggplot(data = expressed_genes_all, aes(x = diet, y = percentage_expressed_genes, fill = location)) +
  geom_bar(stat = "identity", position = 'dodge2') +
  theme(axis.title.y = element_blank()) +
  xlab('Diet') +
  ylab("Percentage of total genes expressed") +
  coord_flip() +
  theme(text = element_text(size=14)) +
  scale_x_discrete(labels = x_labels) +
  theme_light() +
  labs(fill = 'Location') +
  scale_fill_brewer(labels = legend_labels)
dev.off()
