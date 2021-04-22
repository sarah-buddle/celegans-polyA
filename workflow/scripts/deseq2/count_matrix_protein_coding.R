#### Produce count matrix for just protein-coding genes ####

# Load packages
library(rtracklayer)
library(plyranges)

# Import count matrix
# count_matrix <- readRDS('output/deseq2_combined/data/all/count_matrix.rds')
count_matrix <- readRDS(snakemake@input$count_matrix)
count_matrix$gene_id <- rownames(count_matrix)

# Select new and protein coding genes from liftover
liftover_bristol <- rtracklayer::import('from_cluster/liftover_annotations/bristol/liftover_bristol.gtf') %>%
#liftover_bristol <- rtracklayer::import(snakemake@input$liftover_bristol) %>%
  plyranges::filter(type == 'gene' & gene_biotype == 'protein_coding')

liftover_altadena <- rtracklayer::import('from_cluster/liftover_annotations/altadena/liftover_altadena.gtf') %>%
# liftover_altadena <- rtracklayer::import(snakemake@input$liftover_altadena) %>%
  plyranges::filter(type == 'gene' & gene_biotype == 'protein_coding')

bristol_protein_coding <- liftover_bristol$gene_id
altadena_protein_coding <- liftover_altadena$gene_id

protein_coding <- c(bristol_protein_coding, altadena_protein_coding) %>%
  unique(.)

# Only include these genes in count matrix
count_matrix_protein_coding <- plyranges::filter(count_matrix, gene_id %in% protein_coding |
                                                   grepl('WBGene', gene_id) == FALSE) %>%
  plyranges::select(-gene_id)

# Save output
# saveRDS(count_matrix_protein_coding, 'output/deseq2_combined/data/all/count_matrix_protein_coding.rds')
saveRDS(count_matrix_protein_coding, snakemake@output$count_matrix_protein_coding)
