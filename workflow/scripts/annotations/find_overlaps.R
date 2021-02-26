#### Find overlaps between pairs of annotations ####

# Set working directory
setwd("~/OneDrive/Documents/Uni/III/Project/github/celegans-polyA/workflow")

# Load packages
library(GenomicRanges)

# Load granges annotations

#load(file = 'output/annotations/granges/trinity/bristol/as/rep2/trinity_bristol_as_rep2.RData')
#annotation1 <- readRDS(file = 'output/annotations/granges_genes/liftover/altadena/liftover_altadena.rds')
annotation1 <- readRDS(file = snakemake@input$granges1)

#load(file = 'output/annotations/granges/soap/bristol/as/rep2/soap_bristol_as_rep2.rds')
#annotation2 <- readRDS(file = 'output/annotations/granges_genes/soap/altadena/as/rep123/soap_altadena_as_rep123.rds')
annotation2 <- readRDS(file = snakemake@input$granges2)

# Find overlaps
overlaps <- GenomicRanges::findOverlaps(annotation1, annotation2)

# Count overlaps
overlap_counts <- GenomicRanges::countOverlaps(annotation2, annotation1)

# Save outputs
saveRDS(overlaps, file = snakemake@output$overlaps)

saveRDS(overlap_counts, file = snakemake@output$overlap_counts)

# # Testing
# 
# annotation1 <- as.data.frame(sort.GenomicRanges(annotation1))
# head(annotation1, 50)
# annotation2 <- as.data.frame(sort.GenomicRanges(annotation2))
# head(annotation2, 30)
# 
# sum(overlap_counts == 0)
# 
# protein_coding <- subset(annotation1, gene_biotype == 'protein_coding')
# min(protein_coding$width)
# min(annotation2$width)
# 
# plot(protein_coding$width)

# GOI <- subset(data.frame(soap_altadena_as_rep123), ID == 'maker-X-exonerate_est2genome-gene-24.0')
# GOI
# 
# goi_liftover <- subset(data.frame(granges_liftover_altadena), seqnames == 'X' & start == 13557086)
# 
# test <- subset(data.frame(soap_altadena_all_rep123), width > 1000 & width < 10000 & seqnames == 'X')
# test
