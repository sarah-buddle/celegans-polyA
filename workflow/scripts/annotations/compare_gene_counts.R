#### Plot comparison of coverage across annotation types ####

# Load packages
library(rtracklayer)
library(ggplot2)
library(plyranges)

# Load annotations

# Import reference and liftovers
# vc2010 <- readRDS('output/annotations/granges_genes/reference/vc2010/reference_vc2010.rds') %>%
vc2010 <- readRDS(snakemake@input$granges_genes_reference) %>% 
  reduce(.)

# liftover_altadena <- readRDS('output/annotations/granges_genes/liftover/altadena/liftover_altadena.rds') %>%
liftover_altadena <- readRDS(snakemake@input$granges_genes_liftover_altadena) %>% 
  reduce(.)

# liftover_bristol <- readRDS('output/annotations/granges_genes/liftover/bristol/liftover_bristol.rds') %>%
liftover_bristol <- readRDS(snakemake@input$granges_genes_liftover_bristol) %>% 
  reduce(.)

# Make rows fro data frame for reference and liftovers
vc2010_row <- c(nrow(data.frame(vc2010)), 'reference', 'vc2010')
liftover_altadena_row <- c(nrow(data.frame(liftover_altadena)), 'liftover', 'altadena')
liftover_bristol_row <- c(nrow(data.frame(liftover_bristol)), 'liftover', 'bristol')

# Import altadena
# combinedall_altadena <- rtracklayer::import('output/annotations/combined_annotations/combinedall/altadena/combinedall_altadena.gtf') %>%
combinedall_altadena <- rtracklayer::import(snakemake@input$combinedall_altadena) %>% 
  reduce(.)

# combinedmakeronly_altadena <- rtracklayer::import('output/annotations/combined_annotations/combinedmakeronly/altadena/combinedmakeronly_altadena.gtf') %>%
combinedmakeronly_altadena <- rtracklayer::import(snakemake@input$combinedmakeronly_altadena) %>% 
  reduce(.)

# stringtie_altadena <- readRDS('output/annotations/granges/stringtie/altadena/all/rep123/stringtie_altadena_all_rep123.rds') %>%
stringtie_altadena <- readRDS(snakemake@input$granges_genes_stringtie_altadena) %>% 
  plyranges::filter(type == 'transcript') %>%
  reduce(.)

# Make data frame for altadena
gene_counts_altadena <- data.frame(c(nrow(data.frame(combinedall_altadena)),
                            nrow(data.frame(combinedmakeronly_altadena)),
                            nrow(data.frame(stringtie_altadena))))
gene_counts_altadena$annotation_type <- c('all', 'maker', 'stringtie')
colnames(gene_counts_altadena) <- c('genes', 'annotation_type')
gene_counts_altadena$location <- 'altadena'

# Import bristol
# combinedall_bristol <- rtracklayer::import('output/annotations/combined_annotations/combinedall/bristol/combinedall_bristol.gtf')%>%
combinedall_bristol <- rtracklayer::import(snakemake@input$combinedall_bristol) %>% 
    reduce(.)

# combinedmakeronly_bristol <- rtracklayer::import('output/annotations/combined_annotations/combinedmakeronly/bristol/combinedmakeronly_bristol.gtf')%>%
combinedmakeronly_bristol <- rtracklayer::import(snakemake@input$combinedmakeronly_bristol) %>% 
  reduce(.)

# stringtie_bristol <- readRDS('output/annotations/granges/stringtie/bristol/all/rep123/stringtie_bristol_all_rep123.rds') %>%
stringtie_bristol <- readRDS(snakemake@input$granges_genes_stringtie_bristol) %>% 
  plyranges::filter(type == 'transcript')%>%
  reduce(.)

# Make data frame for bristol
gene_counts_bristol <- data.frame(c(nrow(data.frame(combinedall_bristol)),
                                     nrow(data.frame(combinedmakeronly_bristol)),
                                     nrow(data.frame(stringtie_bristol))))
gene_counts_bristol$annotation_type <- c('all', 'maker', 'stringtie')
colnames(gene_counts_bristol) <- c('genes', 'annotation_type')
gene_counts_bristol$location <- 'bristol'


# Combine data frames
gene_counts <- rbind(gene_counts_altadena, gene_counts_bristol) %>%
  rbind(vc2010_row) %>%
  rbind(liftover_altadena_row) %>%
  rbind(liftover_bristol_row)

# Reorder
gene_counts$genes <- as.integer(gene_counts$genes)

gene_counts$annotation_type <- factor(gene_counts$annotation_type,
                                         levels = c('reference', 'liftover',
                                                    'maker', 'stringtie',
                                                    'all'))

gene_counts$location <- factor(gene_counts$location,
                                  levels = c('vc2010', 'bristol', 'altadena'))

# Output table
# write.table(gene_counts, 'output/annotations/gene_counts/gene_counts.tsv')
write.table(gene_counts, snakemake@output$gene_counts)

# Plot
# Labels
x_labels <- c('VC2010', 'N2', 'PS2025')

legend_labels <- c('VC2010 reference', 'Liftover',
                   'MAKER2',
                   'StringTie',
                   'Liftover plus new genes')

# Plot
tiff(snakemake@output$gene_count_plot, width = 1600, height = 800, res = 300)
ggplot2::ggplot(data = gene_counts,
                mapping = aes(x = location, y = genes)) +
  geom_bar(stat = 'identity', aes(fill = annotation_type), position = position_dodge2(width = 0.9, preserve = "single")) +
  scale_fill_brewer(labels = legend_labels,  palette = 'Dark2', direction = -1) +
  theme(text = element_text(size = 10), axis.text.y=element_text(size = 6)) +
  xlab(element_blank()) +
  ylab('Number of unique genes/transcripts') +
  theme_light() +
  theme(panel.grid = element_blank()) +
  theme(legend.key = element_rect(color = NA, fill = NA),
        legend.key.size = unit(0.8, "cm"),
        legend.position="bottom") +
  scale_x_discrete(labels = x_labels) +
  scale_y_continuous(expand = c(0, 0), labels = scales::comma, limits = c(0, 25000)) +
  labs(fill = "Annotation type")
dev.off()

# Save as RDS
p <- last_plot()
saveRDS(p, snakemake@output$gene_count_plot_rds)
