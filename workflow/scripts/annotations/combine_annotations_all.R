# Produce combined annotations (liftver plus new genes) and create list of new genes

# Set working directory
setwd('~/OneDrive/Documents/Uni/III/Project/github/celegans-polyA/workflow')

# Load packages
library(GenomicRanges)
library(plyranges)
library(rtracklayer)

# Import liftover genes
# liftover_genes <- readRDS('output/annotations/granges_genes/liftover/altadena/liftover_altadena.rds')
liftover_genes <- readRDS(snakemake@input$liftover_genes)

# Import new genes from maker and stringtie annotations
# maker_new_genes <- readRDS('output/annotations/combined_annotations/combinedmaker/altadena/new_genes_altadena.rds')
maker_new_genes <- readRDS(snakemake@input$maker_new_genes)

# stringtie_new_genes <- readRDS('output/annotations/combined_annotations/combinedstringtie/altadena/new_genes_altadena.rds')
stringtie_new_genes <- readRDS(snakemake@input$stringtie_new_genes)

# Find new genes common to both annotations
shared_new_genes <- subsetByOverlaps(maker_new_genes, stringtie_new_genes,
                                     ignore.strand = TRUE) %>%
  plyranges::mutate(source = 'maker_stringtie')

#shared_new_genes_df <- data.frame(shared_new_genes)

# New genes only present in one of the annotations
maker_only_new_genes <- IRanges::subsetByOverlaps(maker_new_genes, stringtie_new_genes,
                                         invert = TRUE, ignore.strand = TRUE)

stringtie_only_new_genes <- IRanges::subsetByOverlaps(stringtie_new_genes, maker_new_genes,
                                             invert = TRUE, ignore.strand = TRUE)

all_new_genes <- c(shared_new_genes, stringtie_only_new_genes, maker_only_new_genes) %>%
  plyranges::arrange(seqnames, start) %>%
  plyranges::mutate(type = 'gene')

#saveRDS(all_new_genes, 'output/annotations/combined_annotations/combinedall/altadena/allnewgenes_altadena.rds')
saveRDS(all_new_genes, snakemake@output$all_new_genes)

# all_new_genes_df <- data.frame(all_new_genes)

combined_annotation <- c(liftover_genes, all_new_genes) %>%
  plyranges::arrange(seqnames, start)

# combined_annotation_df <- data.frame(combined_annotation)

# rtracklayer::export(combined_annotation, 'output/annotations/combined_annotations/combinedall/altadena/combinedall_altadena.gtf')
rtracklayer::export(combined_annotation, snakemake@output$all_combined_annotation)
