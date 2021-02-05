#### Make GRanges object maker ####

# Set working directory
setwd("~/OneDrive/Documents/Uni/III/Project/github/celegans-polyA/workflow")

# Load packages
library(rtracklayer)

# Import liftover
# annotation <- rtracklayer::import('from_cluster/maker_annotations/trinity/bristol/as/rep2/trinity_bristol_as_rep2.gff')
annotation <- rtracklayer::import(snakemake@input$maker)

# Save output
save(annotation, file = snakemake@output$granges)