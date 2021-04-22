#### Make GRanges object for liftover annotation ####

# Load packages
library(rtracklayer)
library(plyranges)

# Import annotation
# annotation <- rtracklayer::import('from_cluster/liftover_annotations/altadena/liftover_altadena.gtf')
annotation <- rtracklayer::import(snakemake@input$annotation)

annotation_genes <- subset(annotation, type == 'gene') %>%
  subset(gene_biotype == 'protein_coding') %>%
  unique(.)

# Save output
saveRDS(annotation, file = snakemake@output$granges)
saveRDS(annotation_genes, file = snakemake@output$granges_genes)
