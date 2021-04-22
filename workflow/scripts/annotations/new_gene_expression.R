#### Expression of new genes ####

# Load packages
library(DESeq2)
library(plyranges)

# Load count data
# counts <- data.frame(assay(readRDS('output/deseq2_combined/data/bristol/bristol_full_dds.rds')))
counts <- data.frame(assay(readRDS(snakemake@input$full_dds)))
counts$gene_id <- rownames(counts)
counts <- counts[,c(19, 1:18)]
rownames(counts) <- NULL

# Count data for new genes
new_genes <- counts[grep('(stringtie|maker)', counts$gene_id), ]
options(scipen = 999)
new_genes$mean_counts <- rowMeans(new_genes[2:19])
new_genes <- new_genes[, c(1, 20, 2:19)]

saveRDS(new_genes, snakemake@output$new_gene_counts)

# Load new genes with sources
# allnewgenes_altadena <- readRDS('output/annotations/combined_annotations/combinedall/altadena/allnewgenes_altadena.rds') %>%
#allnewgenes <- readRDS('output/annotations/combined_annotations/combinedall/bristol/allnewgenes_bristol.rds') %>%
allnewgenes <- readRDS(snakemake@input$allnewgenes) %>%
  data.frame(.) %>%
  select(c(gene_id, source))

# Merge to include source column
new_genes <- merge(new_genes, allnewgenes)

# Calculate mean counts for each diet
as_mean_counts <- rowMeans(new_genes[, c(3:5)])
bp_mean_counts <- rowMeans(new_genes[, c(6:8)])
hb101_mean_counts <- rowMeans(new_genes[, c(9:11)])
m9_mean_counts <- rowMeans(new_genes[, c(12:14)])
op50_mean_counts <- rowMeans(new_genes[, c(15:17)])
pf_mean_counts <- rowMeans(new_genes[, c(18:20)])

diet_means <- data.frame(new_genes$gene_id, new_genes$source, as_mean_counts,
                           bp_mean_counts, hb101_mean_counts,m9_mean_counts,
                           op50_mean_counts, pf_mean_counts)

# Count expressed and not expressed genes for a articular annotation type
count_new_genes <- function(annotation_type, sample_means)
{
  all_genes <- nrow(subset(sample_means, new_genes.source == annotation_type))
  expressed_genes <- subset(sample_means, (rowSums(sample_means[3:8] >= 10) > 0)) %>%
    subset(new_genes.source == annotation_type) %>%
    nrow(.)
  not_expressed_genes <- all_genes - expressed_genes
  row1 <- c(annotation_type, 'expressed', expressed_genes)
  row2 <- c(annotation_type, 'not_expressed', not_expressed_genes)
  output <- rbind(row1, row2)
  return(output)
}

# count_new_genes <- function(annotation_type, new_genes_df)
# {
#   all_genes <- nrow(subset(new_genes_df, source == annotation_type))
#   expressed_genes <- nrow(subset(new_genes_df, source == annotation_type &
#                                    mean_counts >= 10))
#                           not_expressed_genes <- all_genes - expressed_genes
#                           row1 <- c(annotation_type, 'expressed', expressed_genes)
#                           row2 <- c(annotation_type, 'not_expressed', not_expressed_genes)
#                           output <- rbind(row1, row2)
#                           return(output)
# }

# Run above function and convert result to data frame
data <- lapply(c('maker', 'StringTie', 'maker_stringtie'), FUN = count_new_genes, sample_means = diet_means) %>%
  do.call(rbind, .) %>% as.data.frame(.)
rownames(data) <- NULL
colnames(data) <- c('annotation_type', 'genes', 'number_of_genes')
data$number_of_genes <- as.integer(data$number_of_genes)

# Export
saveRDS(data, snakemake@output$expression_counts)
