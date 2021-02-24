#### Creates DESeq2 objects for use in further analysis ####

# Set working directory
setwd("~/OneDrive/Documents/Uni/III/Project/github/celegans-polyA/workflow")

# Load packages
library(DESeq2)

# Import sample info
samples <- as.data.frame(read.csv(snakemake@input$samples))

# Only include samples for specified location
samples <- subset(samples, location == snakemake@wildcards$location)

# Make DESeqDataSet object from HTSeq-count data
full_dds <- DESeq2::DESeqDataSetFromHTSeqCount(samples, "from_cluster/htseq_count_export", design = ~ diet)

# Save resulting object
save(full_dds, file = snakemake@output$full_dds)

# Pre-filter data set to remove rows with 0 and 1 counts for all
keep <- rowSums(counts(full_dds)) > 1
dds <- full_dds[keep,]

# Transform data to make sure it is homoskedastic
# rlog better when n < 30, vst better if not
rlog_dds <- DESeq2::rlog(dds, blind = FALSE)

# Save rlog_dds object
save(rlog_dds, file = snakemake@output$rlog_dds)