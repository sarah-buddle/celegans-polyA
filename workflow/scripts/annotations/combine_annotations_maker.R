# Combine Maker annotations for different diets

# Set working directory
setwd('~/OneDrive/Documents/Uni/III/Project/github/celegans-polyA/workflow')

# Load packages
library(GenomicRanges)
library(plyranges)
library(rtracklayer)

# Load data
# liftover <- readRDS('output/annotations/granges/liftover/bristol/liftover_bristol.rds')
liftover <- snakemake@input$granges_liftover

# filepaths <- c('output/annotations/granges_genes/soap/bristol/as/rep123/soap_bristol_as_rep123.rds',
#                'output/annotations/granges_genes/soap/bristol/bp/rep123/soap_bristol_bp_rep123.rds',
#                'output/annotations/granges_genes/soap/bristol/hb101/rep123/soap_bristol_hb101_rep123.rds',
#                'output/annotations/granges_genes/soap/bristol/m9/rep123/soap_bristol_m9_rep123.rds',
#                'output/annotations/granges_genes/soap/bristol/op50/rep123/soap_bristol_op50_rep123.rds',
#                'output/annotations/granges_genes/soap/bristol/pf/rep123/soap_bristol_pf_rep123.rds')
filepaths <- snakemake@input$granges_soap

import_annotations <- function(filepath)
{
  name <- basename(filepath)
  diet <- strsplit(name, '_')[[1]][3]
  annotation <- readRDS(filepath) %>%
    plyranges::mutate(diet = diet) %>%
    plyranges::mutate(midrange = (start + end) / 2 ) %>%
    data.frame(.)
  return(annotation)
}

# op50 <- import_annotations('output/annotations/granges_genes/soap/bristol/op50/rep123/soap_bristol_op50_rep123.rds')
annotations <- lapply(filepaths, FUN = import_annotations)

# finds gene containing midrange position
in_ranges<- function(annotation, this_midrange, chromosome)
{
  matching_range <- subset(annotation, start < this_midrange &
                             end > this_midrange & seqnames == chromosome)
  return(matching_range)
}

# test_in_ranges <- in_ranges(op50, 32077, 'X')

# Initialise combined_annotation data frame
combined_annotation <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(combined_annotation) <- c('chromosome', 'start', 'end', 'diets')

# Produces row for each gene with support from two or more annotations
make_combined_annotation <- function(annotations, this_midrange, chromosome)
{
  matching_genes <- do.call(rbind, lapply(annotations, FUN = in_ranges, this_midrange = this_midrange, chromosome = chromosome))
  if(nrow(matching_genes) > 2) {
    gene <- data.frame(matching_genes$seqnames[1], as.integer(mean(matching_genes$start)),
                       as.integer(mean(matching_genes$end)), I(list(matching_genes$diet)))
    colnames(gene) <- c('chromosome', 'start', 'end', 'diets')
    combined_annotation <<- rbind(combined_annotation, data.frame(gene))
  }
}

# make_combined_annotation(annotations, 32224, 'X')

combine_by_chromosome <- function(annotation_df, chromosome)
{
  chromosome_subset <- subset(annotation_df, seqnames == chromosome)
  lapply(chromosome_subset$midrange, FUN = make_combined_annotation, annotations = annotations,
         chromosome = chromosome)
}

# combine_by_chromosome(op50, 'V')

combine_by_diet <- function(annotation)
{
  chromosomes <- c('I', 'II', 'III', 'IV', 'V', 'X')
  lapply(chromosomes, FUN = combine_by_chromosome, annotation_df = annotation)
}

# combine_by_diet(op50)

# Run nested functions
lapply(annotations, FUN = combine_by_diet)

# Sort and add width column
combined_annotation_full <- arrange(combined_annotation, chromosome, start) %>%
  mutate(width = end - start)

# saveRDS(combined_annotation_full, 'output/annotations/combined_annotations/combinedmaker/bristol/combined_annotation_full_bristol.rds')
saveRDS(combined_annotation_full, snakemake@output$combined_annotation_full)

# Make combined annotation into granges objects
combined_annotation <- GenomicRanges::makeGRangesFromDataFrame(combined_annotation_full,
                                                                  start.field = 'start',
                                                                  end.field='end',
                                                                  ignore.strand = TRUE,
                                                                  keep.extra.columns = TRUE) %>%
  GenomicRanges::reduce(.) %>%
  plyranges::mutate(gene_id = 1:n()) %>%
  plyranges::mutate(gene_id = paste('makergeneid', gene_id, sep = ""))

# rtracklayer::export(combined_annotation, 'output/annotations/combined_annotations/combinedmakeronly/bristol/combinedmakeronly_bristol.gtf')
# saveRDS(combined_annotation, 'output/annotations/combined_annotations/combinedmaker/bristol/combinedmakeronly_bristol.rds')
saveRDS(combined_annotation, snakemake@output$combined_annotation)

# Identify new genes from genes in new combined annotation that don't overlap with anything in the liftover
new_genes <- subsetByOverlaps(combined_annotation, liftover, invert = TRUE, ignore.strand = TRUE) %>%
  plyranges::filter(width < 10000) %>%
  plyranges::mutate(source = 'maker') %>%
  plyranges::mutate(gene_id = 1:n()) %>%
  plyranges::mutate(gene_id = paste('maker', gene_id, sep = ""))

# new_genes_df <- data.frame(new_genes)

# saveRDS(new_genes, 'output/annotations/combined_annotations/combinedmaker/bristol/new_genes_bristol.rds')
saveRDS(new_genes, snakemake@output$new_genes

# Subset genes in liftover annotation
liftover_genes <- filter(liftover, type == 'gene')

# Combine liftover annotation with new genes to make annotation containing the new and existing genes
maker_combined_annotation <- c(liftover_genes, new_genes) %>%
  plyranges::arrange(seqnames, start)

# liftover_combined_annotation_df <- data.frame(liftover_combined_annotation)

# rtracklayer::export(maker_combined_annotation, 'output/annotations/combined_annotations/combinedmaker/bristol/combinedmaker_bristol.gtf')
rtracklayer::export(maker_combined_annotation, snakemake@output$maker_combined_annotation)

# sum(width(liftover_combined_annotation))
# sum(width(liftover_genes))
# sum(width(combined_annotation_gr))
