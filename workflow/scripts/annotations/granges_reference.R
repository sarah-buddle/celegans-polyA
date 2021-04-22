#### Make GRanges object for reference annotation ####


# Load packages
library(rtracklayer)
library(plyranges)

# Import annotation
# annotation <- rtracklayer::import('input/wormbase_annotations/vc2010/vc2010.gff3')
annotation <- rtracklayer::import(snakemake@input$annotation)

annotation_genes <- subset(annotation, type == 'gene') %>%
  subset(biotype == 'protein_coding') %>%
  unique(.)

# Save output
saveRDS(annotation, file = snakemake@output$granges)
saveRDS(annotation_genes, file = snakemake@output$granges_genes)
