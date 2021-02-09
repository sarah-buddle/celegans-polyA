#### Differential expression analysis ####

# set working directory
setwd("~/OneDrive/Documents/Uni/III/Project/github/celegans-polyA/workflow")

# Load packages
library(DESeq2)

# Load dds_full made earlier
load('output/deseq2_data/bristol/bristol_full_dds.RData')
load(snakemake@input$full_dds)

# Perfrom differential expression analysis
ddsDE <- DESeq2::DESeq(full_dds)
dds_res <- DESeq2::results(ddsDE, contrast = c("diet", 
                                                snakemake@wildcards$diet1, 
                                                snakemake@wildcards$diet2))

# Save dds_res object
save(dds_res, file = snakemake@output$dds_res)
