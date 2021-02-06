#### Make GRanges object ####

# Set working directory
setwd("~/OneDrive/Documents/Uni/III/Project/github/celegans-polyA/workflow")

# Load packages
library(rtracklayer)

# Import annotation
# annotation <- rtracklayer::import('from_cluster/maker_annotations/trinity/bristol/as/rep2/trinity_bristol_as_rep2.gff')
# annotation <- rtracklayer::import('from_cluster/liftover_annotations/bristol/liftover_bristol.gtf')
annotation <- rtracklayer::import(snakemake@input$annotation)

# Remove rows corresponding to contigs
annotation <- annotation[which(annotation$type != 'contig'), ]

# Save output
save(annotation, file = snakemake@output$granges)
