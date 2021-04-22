# Produce list of new genes from StringTie annotation

# Load packages
library(GenomicRanges)
library(plyranges)
library(rtracklayer)

# Load annotations
# liftover <- readRDS('output/annotations/granges/liftover/altadena/liftover_altadena.rds')
# stringtie <- readRDS('output/annotations/granges/stringtie/altadena/all/rep123/stringtie_altadena_all_rep123.rds')
liftover <- readRDS(snakemake@input$granges_liftover)
stringtie <- readRDS(snakemake@input$granges_stringtie)

# Find new genes in stringtie annotation
stringtie_new_genes <- IRanges::subsetByOverlaps(stringtie, liftover,
                                                 invert = TRUE,
                                                 ignore.strand = TRUE) %>%
  plyranges::filter(width < 10000) %>%
  plyranges::filter(type == 'transcript') %>%
  plyranges::mutate(gene_id = 1:n()) %>%
  plyranges::mutate(gene_id = paste('stringtie', gene_id, sep = ""))

# stringtie_new_genes_df <- data.frame(stringtie_new_genes)

saveRDS(stringtie_new_genes, snakemake@output$new_genes)

# Add new genes from Stringtie to liftover annotation
liftover_genes <- plyranges::filter(liftover, type == 'gene')

stringtie_combined_annotation <- c(liftover_genes, stringtie_new_genes) %>%
  plyranges::arrange(seqnames, start)

# stringtie_combined_annotation_df <- data.frame(stringtie_combined_annotation)

# rtracklayer::export(stringtie_combined_annotation, 'output/annotations/combined_annotations/combinedstringtie/altadena/combinedstringtie_altadena.gtf')
rtracklayer::export(stringtie_combined_annotation, snakemake@output$stringtie_combined_annotation)
