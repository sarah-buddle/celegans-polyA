#### Differential expression analysis between two diets for same location ####

# Load packages
library(DESeq2)

# Load dds_full made earlier
#load('output/deseq2_data/bristol/bristol_full_dds.RData')
full_dds <- readRDS(snakemake@input$full_dds)

# Perform differential expression analysis
ddsDE <- DESeq2::DESeq(full_dds)
dds_res <- DESeq2::results(ddsDE, contrast = c('diet',
                                                snakemake@wildcards$diet1,
                                                snakemake@wildcards$diet2))

# Save dds_res object
saveRDS(dds_res, file = snakemake@output$dds_res)
