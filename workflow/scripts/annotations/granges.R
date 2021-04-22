#### Make GRanges object ####

# Load packages
library(rtracklayer)

# Import annotation
# annotation <- rtracklayer::import('from_cluster/maker_annotations/trinity/bristol/as/rep2/trinity_bristol_as_rep2.gff')
# annotation <- rtracklayer::import('from_cluster/liftover_annotations/altadena/liftover_altadena.gtf')
annotation <- rtracklayer::import(snakemake@input$annotation)

# Remove rows corresponding to contigs
annotation <- annotation[which(annotation$type != 'contig'), ]

# Just include genes
annotation_genes <- annotation[which(annotation$type == 'gene'), ]

# Save output
saveRDS(annotation, file = snakemake@output$granges)
saveRDS(annotation_genes, file = snakemake@output$granges_genes)
