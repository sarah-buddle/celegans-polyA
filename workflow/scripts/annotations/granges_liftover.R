#### Make GRanges object liftover ####

# Set working directory
setwd("~/OneDrive/Documents/Uni/III/Project/github/celegans-polyA/workflow")

# Load packages
library(rtracklayer)

# Import liftover
liftover <- rtracklayer::import(snakemake@input$liftover)

# Save output
save(liftover, file = snakemake@output$granges)