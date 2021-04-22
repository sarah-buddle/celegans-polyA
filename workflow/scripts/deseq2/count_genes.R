#### Counts unique WB IDs in liftover annotation ####

# Set working directory
setwd("~/OneDrive/Documents/Uni/III/Project/github/celegans-polyA/workflow")

# Load packages
library(rtracklayer)

# Import liftover annotation
# liftover <- rtracklayer::import("from_cluster/liftover_annotations/bristol/liftover_bristol.gtf")
liftover <- rtracklayer::import(snakemake@input$liftover)

# Count unique WB IDs
total_genes <- length(unique(liftover$gene_id))

# Save output
save(total_genes, file = snakemake@output$gene_count)