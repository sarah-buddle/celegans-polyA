#### Calculate coverage of the annotations #### 

# Set working directory
setwd("~/OneDrive/Documents/Uni/III/Project/github/celegans-polyA/workflow")

# Load packages
library(GenomicRanges)

# Load annotation
load(file = snakemake@input$granges)

# Reduce overlapping ranges
reduced_annotation <- GenomicRanges::reduce(annotation)

# Calculate coverage
coverage <- sum(width(reduced_annotation))

# Save output
save(coverage, file = snakemake@output$coverage)