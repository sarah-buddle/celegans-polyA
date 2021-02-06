#### Make GRanges object ####

# Set working directory
setwd("~/OneDrive/Documents/Uni/III/Project/github/celegans-polyA/workflow")

# Load packages
library(rtracklayer)

# Import annotation
annotation <- rtracklayer::import(snakemake@input$annotation)

# Remove rows corresponding to contigs
annotation <- annotation[which(annotation$type != 'contig'), ]

# Save output
save(annotation, file = snakemake@output$granges)