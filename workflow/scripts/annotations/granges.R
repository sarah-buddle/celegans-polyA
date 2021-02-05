#### Make GRanges object ####

# Set working directory
setwd("~/OneDrive/Documents/Uni/III/Project/github/celegans-polyA/workflow")

# Load packages
library(rtracklayer)

# Import annotation
annotation <- rtracklayer::import(snakemake@input$annotation)

# Save output
save(annotation, file = snakemake@output$granges)