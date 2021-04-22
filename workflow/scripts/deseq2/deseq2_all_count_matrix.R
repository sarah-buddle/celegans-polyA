#### Creates count matrix for all data ####

# Set working directory
setwd('~/OneDrive/Documents/Uni/III/Project/github/celegans-polyA/workflow')

# Load packages
library(tidyverse)

# List of paths to samples
# sample_paths <-  c('from_cluster/htseq_count_export/bristol_as_rep1_counts.txt', 
#                    'from_cluster/htseq_count_export/bristol_as_rep2_counts.txt')
sample_paths <- snakemake@input$htseq_count

# imports count data
import_counts <- function(sample_path)
{
  sample_name <- sub('*from_cluster/htseq_count_export/', '', sample_path) %>% 
    sub('_counts.txt', '', .)
  new_column <- readr::read_delim(sample_path, 
                                  col_names = c('gene_id', sample_name), 
                                  delim = '\t')
  return(new_column)
}

# Make matrix of counts
count_matrix <- lapply(sample_paths, import_counts) %>% 
  purrr::reduce(full_join, by = 'gene_id') %>% 
  tibble::column_to_rownames('gene_id')

# Save output
saveRDS(count_matrix, file = snakemake@output$count_matrix)