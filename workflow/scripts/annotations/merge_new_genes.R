#### Combine new gene lists across the diets ####

# Set working directory
setwd('~/OneDrive/Documents/Uni/III/Project/github/celegans-polyA/workflow')

# Load packages
library(GenomicRanges)
library(plyranges)

# Load data file names
# filepaths <- c('output/annotations/granges_new_genes/soap/altadena/as/rep123/soap_altadena_as_rep123.rds', 
#                         'output/annotations/granges_new_genes/soap/altadena/bp/rep123/soap_altadena_bp_rep123.rds')
filepaths <- c(snakemake@input$granges_new_genes)

# granges of overlapping genes between two granges objects
find_common_genes <- function(filepath)
{
  if(exists('common_genes') == TRUE)
  {
    granges_new_genes <- readRDS(filepath)
    overlap_counts <- GenomicRanges::countOverlaps(granges_new_genes, common_genes, 
                                                   minoverlap = 500)
    common_genes <<- plyranges::slice(granges_new_genes, which(overlap_counts != 0))
  }else{
    common_genes <<- readRDS(filepath)
  }
    
}

# find_common_genes('output/annotations/granges_new_genes/soap/altadena/as/rep123/soap_altadena_as_rep123.rds')
# find_common_genes('output/annotations/granges_new_genes/soap/altadena/bp/rep123/soap_altadena_bp_rep123.rds')
# find_common_genes('output/annotations/granges_new_genes/soap/altadena/hb101/rep123/soap_altadena_hb101_rep123.rds')
# find_common_genes('output/annotations/granges_new_genes/soap/altadena/m9/rep123/soap_altadena_m9_rep123.rds')
# find_common_genes('output/annotations/granges_new_genes/soap/altadena/op50/rep123/soap_altadena_op50_rep123.rds')
# find_common_genes('output/annotations/granges_new_genes/soap/altadena/pf/rep123/soap_altadena_pf_rep123.rds')

lapply(filepaths, find_common_genes)

# Save output
saveRDS(common_genes, file = snakemake@output$granges_new_genes)