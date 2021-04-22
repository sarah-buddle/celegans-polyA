#### Comparing numbers of expressed genes between diet treatments ####

# Load packages
library(DESeq2)

# Import sample list
# samples <- as.data.frame(read.csv('scripts/deseq2/htseqcount_samples_full.csv'))
samples <- as.data.frame(read.csv(snakemake@input$samples))

# Set location
samples <-  subset(samples, location == 'altadena')
samples <- subset(samples, location == snakemake@wildcards$location)

# Load appropriate DESeq object made by deseq script
# full_dds <- readRDS('output/deseq2_combined/data/altadena/altadena_full_dds.rds')
full_dds <- readRDS(snakemake@input$full_dds)

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

# count_expressed_genes(full_dds_df, 'altadena_as_rep1', 10)

# apply the function to each sample in our data frame
sample_names <- samples$sample_name
expressed_gene_count <- lapply(sample_names, count_expressed_genes,
                               df = full_dds_df, threshold = 10)

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
saveRDS(expressed_genes, file = snakemake@output$expressed_genes)
