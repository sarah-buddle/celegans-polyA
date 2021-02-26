
# Set working directory
setwd('~/OneDrive/Documents/Uni/III/Project/github/celegans-polyA/workflow')

# Load packages
library(DESeq2)

# Import sample info
#samples <- as.data.frame(read.csv('scripts/deseq2/htseqcount_samples_full.csv'))
samples <- as.data.frame(read.csv(snakemake@input$samples)) 
  
# Only include relevant locations
samples <- subset(samples, location == 'bristol' | location == 'altadena')

# Import count matrix
#count_matrix <- readRDS('output/deseq2_data/all/count_matrix.rds')
count_matrix <- readRDS(snakemake@input$count_matrix)

# Remove rows with extra info rather than gene counts
count_matrix <- subset(count_matrix, grepl('WBGene', rownames(count_matrix)))

# Change NA to 0 or delete rows with NA - doesn't make any different to the PCA plot
count_matrix[is.na(count_matrix)] <- 0
#count_matrix <- na.omit(count_matrix)

# Create deseq object
full_dds <- DESeq2::DESeqDataSetFromMatrix(count_matrix, colData = samples, 
                                           design = ~ diet + location)

# Save resulting object
saveRDS(full_dds, file = snakemake@output$full_dds)

# Transform data to make sure it is homoskedastic
# rlog better when n < 30, vst better if not
vst_dds <- DESeq2::vst(full_dds, blind = FALSE)

# Save rlog_dds object
saveRDS(vst_dds, file = snakemake@output$vst_dds)

# Testing
test <- subset(count_matrix, bristol_as_rep1 == NA)
