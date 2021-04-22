#### Make granges object containing all genes found in new annotation but not liftover ####

# Set working directory
setwd("~/OneDrive/Documents/Uni/III/Project/github/celegans-polyA/workflow")

# Load packages
library(plyranges)

# Load overlap_counts_genes data
#overlap_counts_genes <- readRDS(file = 'output/annotations/overlap_counts_genes/liftover_altadena_vs_soap_altadena_as_rep123.rds')
overlap_counts_genes <- readRDS(file = snakemake@input$overlap_counts_genes)

# Load maker annotation
#granges_genes <- readRDS(file = 'output/annotations/granges_genes/soap/altadena/as/rep123/soap_altadena_as_rep123.rds')
granges_genes <- readRDS(file = snakemake@input$granges_genes)

# Make new granges with just new genes
new_genes <- plyranges::slice(granges_genes, which(overlap_counts_genes == 0))

# Save
saveRDS(new_genes, file = snakemake@output$granges_new_genes)

# Count new genes
new_gene_count <- length(new_genes)

# Write to table

#sample_name <- 'soap_bristol_as_rep2'
sample_name <- paste(snakemake@wildcards$annotation_type, '_',
                     snakemake@wildcards$location, '_', 
                     snakemake@wildcards$diet, '_', 
                     snakemake@wildcards$replicate, sep = "")

# New row
# new_row <- c(sample_name, 'soap', 'bristol', 'as', 'rep2', 100)
new_row <- c(sample_name, snakemake@wildcards$annotation_type, 
             snakemake@wildcards$location, snakemake@wildcards$diet, 
             snakemake@wildcards$replicate, new_gene_count)

# Making or updating coverage_table.txt
# Checks whether coverage_table.txt file already exists
if(file.exists('output/annotations/new_genes/new_gene_counts.txt') == TRUE) {
  
  # import existing table
  new_gene_counts <- read.table('output/annotations/new_genes/new_gene_counts.txt',
                               header = TRUE)
  
  # removes row corresponding to current sample if necessary
  new_gene_counts <- new_gene_counts[new_gene_counts$sample_name != sample_name, ]
  
  # adds new row
  new_gene_counts <- rbind(new_gene_counts, new_row)
  
  # sort alphabetically by sample name
  new_gene_counts <- new_gene_counts[order(new_gene_counts$sample_name), ]
  
  # write  to new_gene_counts.txt
  write.table(new_gene_counts, file = 'output/annotations/new_genes/new_gene_counts.txt',
              row.names = FALSE)
  
  # if no coverage_table.txt file 
}else{
  # column headings
  headers <- c('sample_name', 'annotation_type', 'location', 'diet', 'replicate', 'new_gene_count')
  
  # Makes data frame from headings nd new row
  new_gene_counts <- t(data.frame(headers, new_row))
  
  # Write to coverage_table.txt
  write.table(new_gene_counts, file = 'output/annotations/new_genes/new_gene_counts.txt', 
              col.names = FALSE, row.names = FALSE)
}