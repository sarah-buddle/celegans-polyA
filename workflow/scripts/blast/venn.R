# Venn diagram to display homology analysis

# Set working directory
setwd("~/OneDrive/Documents/Uni/III/Project/github/celegans-polyA/workflow")

# Install ggvenn (not available on conda)
install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")

# Load packages
library(tidyverse)
library(ggvenn)
library(ggpubr)

# Load data
# new_gene_counts <- readRDS('output/annotations/new_gene_expression/data/bristol/new_gene_counts_bristol.rds')
new_gene_counts <- readRDS(snakemake@input$new_gene_counts)

# no_hits_genome <- read.table('output/blast/vc2010_genome_no_hits/bristol/names_vc2010_genome_no_hits_bristol.txt')  %>% 
no_hits_genome <- read.table(snakemake@input$no_hits_genome) %>% 
  dplyr::rename(gene_id = V1)

# no_hits_transcripts <- read.table('output/blast/vc2010_transcripts_no_hits/bristol/names_vc2010_transcripts_no_hits_bristol.txt') %>% 
no_hits_transcripts <- read.table(snakemake@input$no_hits_transcripts) %>% 
  dplyr::rename(gene_id = V1)

# no_hits_proteins <- read.table('output/blast/vc2010_proteins_no_hits/bristol/names_vc2010_proteins_no_hits_bristol.txt') %>% 
no_hits_proteins <- read.table(snakemake@input$no_hits_proteins) %>% 
  dplyr::rename(gene_id = V1)

# hits_related_species_transcripts <- read.table('output/blast/related_species_transcripts/bristol/names_related_species_transcripts_hits_bristol.txt') %>% 
hits_related_species <- read.table(snakemake@input$hits_related_species) %>% 
  dplyr::rename(gene_id = V1)

# hits_related_species_proteins <- read.table('output/blast/related_species_proteins/bristol/names_related_species_proteins_hits_bristol.txt') %>% 
hits_related_species <- read.table(snakemake@input$hits_related_species) %>% 
  dplyr::rename(gene_id = V1)

# Calculate diet means and count genes as expressed if at least 1 has mean read count of greater than 10
mean_cols <- function(x, col1, col2, col3){
  (x[[col1]] + x[[col2]] + x[[col3]]) / 3
}

new_gene_counts_adj <- new_gene_counts %>% 
  mutate(mean_as = mean_cols(., 3, 4, 5)) %>% 
  mutate(mean_bp = mean_cols(., 6, 7, 8)) %>% 
  mutate(mean_hb101 = mean_cols(., 9, 10, 11)) %>% 
  mutate(mean_m9 = mean_cols(., 12, 13, 14)) %>% 
  mutate(mean_op50 = mean_cols(., 15, 16, 17)) %>% 
  mutate(mean_pf = mean_cols(., 18, 19, 20)) %>% 
  mutate(expressed = if_else(mean_as > 10 | mean_bp > 10 | mean_hb101 > 10 |
                               mean_m9 > 10 | mean_op50 > 10 | mean_pf > 10, TRUE, FALSE)) %>% 
  select(c(gene_id, expressed)) %>% 
  mutate(gene_id = paste('bristol', gene_id, sep = '_')) %>% 
  mutate(no_hits_genome = gene_id %in% no_hits_genome$gene_id) %>% 
  mutate(no_hits_transcripts = gene_id %in% no_hits_transcripts$gene_id) %>% 
  mutate(no_hits_proteins = gene_id %in% no_hits_proteins$gene_id) %>% 
  mutate(hits_related_species_transcripts = gene_id %in% hits_related_species_transcripts$gene_id) %>% 
  mutate(hits_related_species_proteins = gene_id %in% hits_related_species_proteins$gene_id)

# Make table

type <- rep(c('all_new_genes', 'no_hits_genome', 'no_hits_transcripts', 
           'no_hits_proteins', 'hits_related_species_transcripts', 
           'hits_related_species_proteins'), each = 2)

expressed <- rep(c('expressed', 'low_expression'), times = 6)

gene_counts <- c(nrow(subset(new_gene_counts_adj, expressed == TRUE)),
                 nrow(subset(new_gene_counts_adj, expressed == FALSE)),
                 nrow(subset(new_gene_counts_adj, expressed == TRUE & no_hits_genome == TRUE )),
                 nrow(subset(new_gene_counts_adj, expressed == FALSE & no_hits_genome == TRUE )),
                 nrow(subset(new_gene_counts_adj, expressed == TRUE & no_hits_transcripts == TRUE )),
                 nrow(subset(new_gene_counts_adj, expressed == FALSE & no_hits_transcripts == TRUE )),
                 nrow(subset(new_gene_counts_adj, expressed == TRUE & no_hits_proteins == TRUE )),
                 nrow(subset(new_gene_counts_adj, expressed == FALSE & no_hits_proteins == TRUE )),
                 nrow(subset(new_gene_counts_adj, expressed == TRUE & hits_related_species_transcripts == TRUE )),
                 nrow(subset(new_gene_counts_adj, expressed == FALSE & hits_related_species_transcripts == TRUE )),
                 nrow(subset(new_gene_counts_adj, expressed == TRUE & hits_related_species_proteins == TRUE )),
                 nrow(subset(new_gene_counts_adj, expressed == FALSE & hits_related_species_proteins == TRUE )))

count_table <- data.frame(type, expressed, gene_counts)
saveRDS(count_table, 'output/blast/count_table/bristol/count_table_bristol.rds')

# Make list for Venn diagram
venn_list <- list(
  expressed = subset(new_gene_counts_adj, expressed == TRUE) %>% 
    .[['gene_id']],
  no_hits_genome = subset(new_gene_counts_adj, no_hits_genome == TRUE) %>%
    .[['gene_id']],
  no_hits_transcripts = subset(new_gene_counts_adj, no_hits_transcripts == TRUE) %>% 
    .[['gene_id']],
  no_hits_proteins = subset(new_gene_counts_adj, no_hits_proteins == TRUE) %>% 
    .[['gene_id']]
  # hits_related_species = subset(new_gene_counts_adj, hits_related_species == TRUE) %>% 
  #   .[['gene_id']]
)

# Plot venn diagram
# tiff('output/blast/venn/bristol/venn_bristol.tiff', width = 1600, height = 800, res = 300)
tiff(snakemake@output$venn, width = 1600, height = 800, res = 300)
ggvenn(venn_list, stroke_size = 0.1, set_name_size = 2, text_size = 1.5, 
       fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"))
dev.off()


# Plot legend
x <- c(1:4)
y <- x
z <- c('Expressed', 'No hits genome',  'No hits transcripts', 'No hits proteins')

df <- data.frame(x, y, z)

df$z <- factor(df$z, levels = c('Expressed', 'No hits genome',  'No hits transcripts', 'No hits proteins'))

legend_labels <- c('\nHigher expression          \n', '\nNo hits to\nVC2010 genome          \n',
                   '\nNo hits to\nVC2010 transcripts          \n', '\nNo hits to\nVC2010 proteins          \n')

p <- ggplot(data = df, aes(x, y, fill = z)) +
  geom_bar(stat = 'identity', alpha = 0.5) +
  scale_fill_manual(values = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"), 
                    labels = legend_labels) +
  theme(legend.position = 'bottom', legend.title = element_blank(), 
        text = element_text(size = 16), legend.text.align = 0.3) +
  guides(fill=guide_legend(ncol=2))

tiff(snakemake@output$legend, width = 1600, height = 800, res = 300)
ggpubr::get_legend(p) %>% 
  as_ggplot(.)
dev.off()




