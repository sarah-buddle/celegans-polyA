#### Creates DESeq2 object without m9/starvation diet treatment ####

# Load packages
library(DESeq2)

# Import sample info
samples <- as.data.frame(read.csv(snakemake@input$samples))

# Only include samples from a specifed location
samples <- subset(samples, location == snakemake@wildcards$location)

# Remove starvation
samples <- subset(samples, diet != "m9")

# Make DESeqDataSet object from HTSeq-count data
full_dds <- DESeq2::DESeqDataSetFromHTSeqCount(samples, "from_cluster/htseq_count_export",
                                               design = ~ diet)

# pre-filter data set to remove rows with 0 and 1 counts for all
keep <- rowSums(counts(full_dds)) > 1
dds <- full_dds[keep,]

# transform data to make sure it is homoskedastic
# rlog better when n < 30, vst better if not
rlog_dds <- DESeq2::rlog(dds, blind = FALSE)

# save rlog_dds object
saveRDS(rlog_dds, file = snakemake@output$rlog_dds)
