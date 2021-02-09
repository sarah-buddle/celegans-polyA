#### Find overlaps between pairs of annotations ####

# Set working directory
setwd("~/OneDrive/Documents/Uni/III/Project/github/celegans-polyA/workflow")

# Load packages
library(GenomicRanges)

# Load granges annotations

load(file = snakemake@input$granges1)
annotation1 <- annotation

load(file = snakemake@input$granges2)
annotation2 <- annotation

# Find overlaps
overlaps <- GenomicRanges::findOverlaps(annotation1, annotation2)

# Save output
save(overlaps, file = snakemake@output$overlaps)
