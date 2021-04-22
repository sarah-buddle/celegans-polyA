# Plot of new genes with homology to transcripts and proteins in related species

# Load packages
library(tidyverse)

# Load count tables
# counts_bristol <- readRDS('output/blast/count_table/bristol/count_table_bristol.rds') %>%
counts_bristol <- readRDS(snakemake@input$counts_bristol) %>%
  dplyr::mutate(location = 'bristol')
# counts_altadena <- readRDS('output/blast/count_table/altadena/count_table_altadena.rds') %>%
counts_altadena <- readRDS(snakemake@input$counts_altadena) %>%
  dplyr::mutate(location = 'altadena')

# Combine
counts <- rbind(counts_bristol, counts_altadena) %>%
  subset(type == 'hits_related_species_proteins' | type == 'hits_related_species_transcripts')

# Reorder
counts$type <- factor(counts$type, levels = c('hits_related_species_transcripts',
                                              'hits_related_species_proteins'))

counts$expressed <- factor(counts$expressed, levels = c('low_expression', 'expressed'))

counts$location <- factor(counts$location, levels = c('bristol', 'altadena'))

# Labels
x_labels <- c('Transcripts', 'Proteins')
legend_labels <- c('Very low expression', 'Higher expression')

locations <- c('N2', 'PS2025')
names(locations) <- c('bristol', 'altadena')

# Plot
# tiff('output/blast/related_species_plot/related_species_plot.tiff', width = 1600, height = 800, res = 300)
tiff(snakemake@output$related_species_plot, width = 1600, height = 800, res = 300)
ggplot2::ggplot(data = counts, mapping = aes(x = type, y = gene_counts, fill = expressed)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~ location, labeller = labeller(location = locations)) +
  scale_fill_brewer(labels = legend_labels) +
  ylab('Number of predicted genes') +
  theme_light() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key = element_rect(color = NA, fill = NA),
        legend.key.size = unit(0.5, "cm"),
        legend.position = c(0.14, 0.84),
        axis.title.x = element_blank()) +
  scale_x_discrete(labels = x_labels) +
  scale_y_continuous(expand = c(0, 0), labels = scales::comma, limits = c(0, 40))
dev.off()
