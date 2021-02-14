#### Comparing numbers of expressed genes between diet treatments ####

# Set working directory
setwd("~/OneDrive/Documents/Uni/III/Project/github/celegans-polyA/workflow")

# Load packages
library(DESeq2)

# Import sample list
samples <- as.data.frame(read.csv(snakemake@input$samples))

# Set location
samples <- subset(samples, location == snakemake@wildcards$location)

# Load appropriate DESeq object made by deseq script
# load('output/deseq2_data/altadena/altadena_full_dds.RData')
load(snakemake@input$full_dds)

# convert DESeq object into data frame
full_dds_df <- as.data.frame(assay(full_dds))
full_dds_df$gene_id <- rownames(full_dds_df)
rownames(full_dds_df) <- NULL

# Count total unique gene IDs
total_genes <- length(unique(full_dds_df$gene_id))

# function to count genes with reads above a certain threshold
count_expressed_genes <- function(df, sample, threshold)
{
  output <- sum(df[, sample] > threshold)
  return(output)
}

# apply the function to each sample in our data frame
sample_names <- samples$sample_name
expressed_gene_count <- lapply(sample_names, count_expressed_genes, 
                               df = full_dds_df, threshold = 1)

# make data frame of results
expressed_genes <- do.call(rbind, Map(data.frame, sample_name = sample_names, 
                                      expressed_gene_count=expressed_gene_count,
                                      location=samples$location,
                                      diet= samples$diet))
rownames(expressed_genes) <- NULL

# Calculate percentage of total genes expressed
expressed_genes$percentage_expressed_genes <- 
  expressed_genes$expressed_gene_count*100/total_genes

# Save result
save(expressed_genes, file = snakemake@output$expressed_genes)
