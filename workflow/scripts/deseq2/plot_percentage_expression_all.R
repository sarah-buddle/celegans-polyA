#### Plot percentage expressed genes across all samples ####

# Load packages
library(ggplot2)

# Load data objects
# expressed_genes1 <- readRDS('output/deseq2_combined/expressed_genes/altadena/altadena_expressed_genes.rds')
expressed_genes1 <- readRDS(snakemake@input$full_dds_altadena)

# expressed_genes2 <- readRDS('output/deseq2_combined/expressed_genes/bristol/bristol_expressed_genes.rds')
expressed_genes2 <- readRDS(snakemake@input$full_dds_bristol)

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
                               levels = c("as", "bp", "hb101", "op50", "pf", "m9"))
expressed_genes_all$location <- factor(expressed_genes_all$location,
                                       levels = c('bristol', 'altadena'))



x_labels <- c('AS', 'BP', 'HB101', 'OP50', 'PF', 'Starvation')

locations <- c('N2', 'PS2025')
names(locations) <- c('bristol', 'altadena')

# Plot
#tiff('output/deseq2_plots/all/percentage_plot_all.tiff', width = 1600, height = 800, res = 300)
tiff(snakemake@output$percentage_plot, width = 1600, height = 800, res = 300)
ggplot(data = expressed_genes_all, aes(x = diet, y = percentage_expressed_genes, fill = diet)) +
  geom_bar(stat = "identity", position = 'dodge2') +
  facet_wrap(~ location, labeller = labeller(location = locations)) +
  xlab('Diet') +
  ylab("Percentage of total genes expressed") +
  scale_x_discrete(labels = x_labels) +
  theme_light() +
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(size = 6)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 40)) +
  scale_fill_brewer(palette = 'Set2')
  guides(fill = guide_legend(reverse = TRUE))
dev.off()
