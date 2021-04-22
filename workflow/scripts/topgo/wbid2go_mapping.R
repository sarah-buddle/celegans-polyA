#### Create file mapping WB IDs to GO terms ####

# set working directory
setwd("~/OneDrive/Documents/Uni/III/Project/github/celegans-polyA/workflow")

# Load packages
library(tidyverse)

# import mapping of gene names to GO terms
GO <- readr::read_tsv(snakemake@input$go_terms, col_names = FALSE) %>% 
  .[c("X2", "X5")] %>% # keep just gene id and associated GO term
  unique(.) # get rid of duplicates
colnames(GO) <- c("gene_id", "GO_term") # rename columns

# function to create list of GO terms in format required for topGO
collapse_GO <- function(gene)
{
  a <- subset(GO, gene_id == gene)
  output <- paste(a$GO_term, collapse=", ")
}

# all the gene_ids that we have GO terms for
gene_ids <- unique(GO$gene_id)

# use the collapse_GO function
geneID2GO_prep <- sapply(gene_ids, collapse_GO) %>% 
  as.data.frame(.)

# write to file that can be read by readMappings
write.table(geneID2GO_prep, file = snakemake@output$wbid2go, row.names = TRUE,
            col.names = FALSE, quote = FALSE, sep  = '\t')
