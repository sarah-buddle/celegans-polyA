# Set working directory
setwd('~/OneDrive/Documents/Uni/III/Project/github/celegans-polyA/workflow')

# Load packages
library(tidyverse)

# Load count matrix for all
count_matrix <- readRDS('output/deseq2_combined/data/all/count_matrix_protein_coding.rds')
count_matrix <- readRDS(snakemake@input$count_matrix)

# Get rid of unaligned reads
count_matrix <- subset(count_matrix, grepl('WBGene', rownames(count_matrix)) |
                         grepl('maker', rownames(count_matrix)) |
                         grepl('stringtie', rownames(count_matrix)))

#Load file paths
filepaths <- c('output/deseq2_combined/results_location/all/all_bristol_altadena.rds',
               'output/deseq2_combined/results/bristol/bristol_m9_op50.rds',
               'output/deseq2_combined/results/bristol/bristol_m9_pf.rds',
               'output/deseq2_combined/results/bristol/bristol_pf_op50.rds',
               'output/deseq2_combined/results/bristol/bristol_pf_as.rds',
               'output/deseq2_combined/results/bristol/bristol_as_op50.rds',
               'output/deseq2_combined/results/altadena/altadena_m9_op50.rds',
               'output/deseq2_combined/results/altadena/altadena_m9_pf.rds',
               'output/deseq2_combined/results/altadena/altadena_pf_op50.rds',
               'output/deseq2_combined/results/altadena/altadena_pf_as.rds',
               'output/deseq2_combined/results/altadena/altadena_as_op50.rds'
)
filepaths <- snakemake@input$dds_res

# Identify genes present in one but not both isolates
accessory_genes_all <- count_matrix[!complete.cases(count_matrix), ]

accessory_genes <- rownames(accessory_genes_all)

accessory_genes_bristol <- accessory_genes_all %>% 
  subset(is.na(altadena_as_rep1)) %>% 
  rownames(.)

accessory_genes_altadena <- accessory_genes_all %>% 
  subset(is.na(bristol_as_rep1)) %>% 
  rownames(.)

# all proportion of accessory genes
all <- length(accessory_genes) / nrow(count_matrix)

bristol_all <- length(accessory_genes_bristol) / 
  nrow(subset(count_matrix, is.na(bristol_as_rep1) == FALSE))

altadena_all <- length(accessory_genes_altadena) / 
  nrow(subset(count_matrix, is.na(altadena_as_rep1) == FALSE))

# Make data frame

proportions <- data.frame(c('all', 'bristol_all', 'altadena_all'), 
                          c(all, bristol_all, altadena_all))
colnames(proportions) <- c('name', 'proportion_accessory')

proportion_accessory_genes <- function(filepath)
{
  # comparison
  name <- basename(filepath) %>% 
    stringr::str_remove('.rds')
  
  # Load dds_res object
  dds_res <- data.frame(readRDS(filepath))
  
  # Subset differentially expressed genes
  differentially_expressed_genes <- rownames(subset(dds_res, padj < 0.05))
  
  # Calculate proportion of accessory genes in differentially expressed genes
  DE <- sum(accessory_genes %in% differentially_expressed_genes) / 
    length(differentially_expressed_genes)
  
  # Append to proportions data frame
  new_row <- c(name, DE)
  proportions <<- rbind(proportions, new_row)
}

lapply(filepaths, proportion_accessory_genes)
proportions$proportion_accessory <- as.numeric(proportions$proportion_accessory)

proportions <- cbind(proportions, str_split_fixed(proportions$name, '_', 2))
names(proportions)[names(proportions) == '1'] <- 'location'
names(proportions)[names(proportions) == '2'] <- 'comparison'
proportions$comparison[proportions$name == 'all'] <- 'all'


# Reorder

proportions$comparison <- factor(proportions$comparison,
                                 levels = c('all', 'bristol_altadena', 
                                            'm9_op50', 'm9_pf', 'pf_as', 
                                            'pf_op50', 'as_op50'))

proportions$location <- factor(proportions$location, 
                               levels = c('all', 'bristol', 'altadena'))

proportions_total <- subset(proportions, comparison == 'all')
proportions_de <- subset(proportions, comparison != 'all')

# Labels and colors

legend_labels <- c("All", "N2", "PS2025")

#library(viridis)
#colors <- viridis(5, option = 'viridis')[2:4]
colors <- c("#3B528BFF", "#21908CFF", "#5DC863FF")

x_labels <- c('N2 -\nPS2025', 'Starvation\n- OP50', 'Starvation\n- PF', 
              'PF - OP50', 'PF - AS', 'AS - OP50')

# Plot

#tiff(snakemake@output$proportion_new_genes_plot, width = 1600, height = 800, res = 300)
tiff('output/deseq2_combined/plots/proportion_new_genes/3/proportion_new_genes.tiff', width = 1600, height = 800, res = 300)
ggplot2::ggplot(data = proportions_de) +
  geom_bar(stat = 'identity', position = position_dodge2(width = 0.9, preserve = "single"), 
           mapping = aes(x = comparison, y = 100*proportion_accessory, 
                         fill = location)) +
  geom_hline(proportions_total,
             mapping = aes(yintercept = 100*proportion_accessory, color = location), 
             linetype = 'dashed', size = 0.5) +
  annotate('text', x = 6.35, y = 9.6, label = 'All', size = 3) +
  annotate('text', x = 6.3, y = 6.1, label = 'N2', size = 3) +
  annotate('text', x = 6.1, y = 4.45, label = 'PS2025', size = 3) +
  xlab(NULL) +
  ylab('Percentage accessory genes') +
  scale_x_discrete(labels = x_labels) +
  theme_light() +
  scale_color_manual(labels = legend_labels, values = colors) +
  scale_fill_manual(labels = legend_labels, values = colors) +
  labs(fill = 'Strain') +
  theme(panel.grid = element_blank(), text = element_text(size = 9), 
        legend.position = 'right') +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 12)) +
  guides(color = FALSE)
dev.off()


# Accessory genes plot
altadena_total <- nrow(subset(count_matrix, is.na(altadena_as_rep1) == FALSE))
accessory_genes_altadena_liftover <- sum(grepl('WBGene', accessory_genes_altadena))
accessory_genes_altadena_transcriptome <- sum(grepl('maker', accessory_genes_altadena)) +
  sum(grepl('ring', accessory_genes_altadena))

bristol_total <- nrow(subset(count_matrix, is.na(bristol_as_rep1) == FALSE))
accessory_genes_bristol_liftover <- sum(grepl('WBGene', accessory_genes_bristol))
accessory_genes_bristol_transcriptome <- sum(grepl('maker', accessory_genes_bristol)) +
  sum(grepl('ring', accessory_genes_altadena))

total <- c(bristol_total, altadena_total)
liftover_accessory <- c(accessory_genes_bristol_liftover, accessory_genes_altadena_liftover)
transcriptome_accessory <- c(accessory_genes_bristol_transcriptome, accessory_genes_altadena_transcriptome)

location <- c('bristol', 'altadena')

accessory_plot <- data.frame(location, total, liftover_accessory, transcriptome_accessory) %>%
  mutate(shared = total - liftover_accessory - transcriptome_accessory) %>%
  select(-total) %>%
  pivot_longer(!location, names_to = 'type', values_to = 'number_of_genes')

# Reorder
accessory_plot$location <- factor(accessory_plot$location, levels = c('bristol', 'altadena'))
accessory_plot$type <- factor(accessory_plot$type, levels = c('transcriptome_accessory',
'liftover_accessory', 'shared'))

# Labels
x_labels <- c('N2', 'PS2025')
legend_labels <- c('Accessory genes from\ntranscriptome-based annotations',
                   'Accessory genes from\nliftover annotation',
                   'Shared genes')

tiff(snakemake@output$accessory_genes_plot, width = 1600, height = 800, res = 300)
ggplot2::ggplot(data = accessory_plot, mapping = aes(x = location, y = number_of_genes,
                                                     fill = type)) +
  geom_bar(stat = 'identity') +
  theme_light() +
  theme(panel.grid = element_blank(),
        legend.key.size = unit(1, "cm"),
        text = element_text(size = 10)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 25000)) +
  scale_fill_brewer(palette = 'Set2', labels = legend_labels) +
  xlab('Strain') +
  ylab('Number of genes') +
  labs(fill = 'Type of gene') +
  scale_x_discrete(labels = x_labels)
dev.off()






