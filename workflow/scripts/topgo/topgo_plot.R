#### Plot of GO analysis ####

# Load packages
library(tidyverse)

# Load data
# master_gentable <- readRDS('output/topgo/master_gentable/1/master_gentable.rds')
master_gentable <- readRDS(snakemake@input$master_gentable)
# rrvgo <- readRDS('output/topgo/rvvgo/rrvgo.rds')
rrvgo <- readRDS(snakemake@input$rrvgo)

# Identify GO terms with highest mean padj
top_go_terms <- aggregate(master_gentable, by = list(master_gentable$GO.ID), FUN = mean) %>%
  top_n(100, wt = logpadj)

# Only look at these GO terms
master_gentable <- subset(master_gentable, GO.ID %in% top_go_terms$Group.1)

# Merge master_gentable and rrvgo, sort by parentTerm
full_data <- merge.data.frame(master_gentable, rrvgo, by = 'GO.ID') %>%
  dplyr::arrange(parentTerm)

# Reorder
full_data$parentTerm <- factor(full_data$parentTerm, levels = rev(unique(as.character(full_data$parentTerm))))
full_data$term <- factor(full_data$term, levels = rev(unique(as.character(full_data$term))))
full_data$name <- factor(full_data$name, levels = c('bristol_m9_op50', 'altadena_m9_op50',
                                                    'bristol_m9_pf', 'altadena_m9_pf',
                                                    'bristol_pf_op50', 'altadena_pf_op50',
                                                    'bristol_pf_as', 'altadena_pf_as',
                                                    'bristol_as_op50', 'altadena_as_op50'))
full_data$comparison <- factor(full_data$comparison, levels = c('m9_op50',
                                                          'm9_pf','pf_op50',
                                                          'pf_as','as_op50'))

# Labels

x_labels <- rep(c('N2', 'PS2025'), times = 5)

facet_labels <- c('Starvation\n - OP50', 'Starvation\n - PF', 'PF - OP50',
                  'PF - AS', 'AS - OP50')
names(facet_labels) <- c('m9_op50', 'm9_pf', 'pf_op50', 'pf_as','as_op50')

# Plot
tiff(snakemake@output$topgo_plot, width = 2200, height = 2500, res = 300)
# tiff('output/topgo/plots/1/topgo_plot_1.tiff', width = 2200, height = 2500, res = 300)
ggplot2::ggplot(full_data, mapping = aes(x = location, y = term, fill = logpadj)) +
  geom_tile() +
  theme_light() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_gradient(guide=guide_colourbar(reverse = FALSE), low='white', high='navy') +
  labs(fill = 'Enrichment') +
  ylab('GO term') +
  theme(axis.title.x = element_blank(), axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 7)) +
  scale_x_discrete(labels = x_labels) +
  facet_grid(. ~ comparison, labeller = labeller(comparison = facet_labels))
dev.off()
