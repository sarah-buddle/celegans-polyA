# from Stringtie manual and vignettes, both in Zotero

#### Setup ####

# install packages
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()

#BiocManager::install("DESeq2") # don't update rlang at this point
#install.packages("gplots")
#install.packages("RColorBrewer")
#install.packages("tidyverse")
#install.packages("dplyr")

# import libraries
library("DESeq2")
library( "gplots" )
library( "RColorBrewer" )
# library("genefilter")
library("tidyverse")
library("dplyr")

# set working directory
setwd("~/OneDrive/Documents/Uni/III/Project/Scripts/polyA/genome_mapping")

#### Import and Check Data ####

# load gene and transcript count matrices
geneCountData <- as.matrix(read.csv("gene_count_matrix.csv", row.names="gene_id"))
transcriptCountData <- as.matrix(read.csv("transcript_count_matrix.csv", row.names="transcript_id"))

# load labels
colData <- read.csv("samples.csv", row.names=1)

# check all sample IDs in colData are also in CountData and match their orders
# genes
all(rownames(colData) %in% colnames(geneCountData))
geneCountData <- geneCountData[, rownames(colData)]
all(rownames(colData) == colnames(geneCountData))

# transcripts
all(rownames(colData) %in% colnames(transcriptCountData))
transcriptCountData <- transcriptCountData[, rownames(colData)]
all(rownames(colData) == colnames(transcriptCountData))

#### Make DESeq Object ####

# create a DESeqDataSet from count matrix and labels and make DESeq object
# genes
dds_genes <- DESeqDataSetFromMatrix(countData = geneCountData,
                              colData = colData, design = ~ diet) %>% 
  DESeq(.)

# transcripts
dds_transcripts <- DESeqDataSetFromMatrix(countData = transcriptCountData,
                                    colData = colData, design = ~ diet) %>% 
  DESeq(.)

#### General Plots ####

# genes

# log transform
rlog_dds_genes <- rlog(dds_genes)

# PCA plot
plotPCA(rlog_dds_genes, intgroup = "diet")
PCA_genes <- plotPCA(rlog_dds_genes, intgroup = "diet", returnData = TRUE)

# heat map
sampleDistGenes <- dist( t( assay(rlog_dds_genes) ) )
sampleDistMatrixGenes <- as.matrix( sampleDistGenes )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrixGenes, trace="none", col=colours)

# gene clustering
topVarGenes <- head( order( rowVars( assay(rlog_dds_genes) ), decreasing=TRUE ), 35 )
heatmap.2( assay(rlog_dds_genes)[ topVarGenes, ], scale="row",
           trace="none", dendrogram="column",
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))


# transcripts

# log transform
rlog_dds_transcripts <- rlog(dds_transcripts)

# PCA plot
plotPCA(rlog_dds_transcripts, intgroup = "diet")
PCA_transcripts <- plotPCA(rlog_dds_transcripts, intgroup = "diet", returnData = TRUE)

# heat map
sampleDistTranscripts <- dist( t( assay(rlog_dds_transcripts) ) )
sampleDistMatrixTranscripts <- as.matrix( sampleDistTranscripts )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrixTranscripts, trace="none", col=colours)

# transcript clustering
topVarTranscripts <- head( order( rowVars( assay(rlog_dds_transcripts) ), decreasing=TRUE ), 35 )
heatmap.2( assay(rlog_dds_transcripts)[ topVarTranscripts, ], scale="row",
           trace="none", dendrogram="column",
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))


#### Analysis of Specific Pairs ####

# m9 vs op50

# generate results tables
res_genes_m9vsop50 <- results(dds_genes, contrast = c("diet", "m9", "op50"))

plotMA(res_genes_m9vsop50)

# subset genes/transcripts with p < 0.1
resSigGenesm9vsop50 <- res_genes_m9vsop50[which(res_genes_m9vsop50$padj < 0.1), ]
nrow(geneCountData) # total number of genes
nrow(resSigGenesm9vsop50) # number of genes with p < 0.1

# export for gene ontology analysis

# function to cut off last part of gene ID, so it can be used by the online tool
WBgeneID <- function(geneID)
{ return(substr(geneID, start = 1, stop = 14))
}

# take just geneIDs from significant transcripts above and apply above function
resSigNamesGenesm9vsop50 <- as.data.frame(rownames(resSigGenesm9vsop50)) %>% 
  lapply(., WBgeneID)

# write to csv for use in online tool
# https://wormbase.org/tools/enrichment/tea/tea.cgi
write.csv(resSigNamesGenesm9vsop50, file="m9vsop50sig.csv", row.names = FALSE)


# find genes with largest fold change with padj < 0.1
# strongest downregulation
head( resSigGenesm9vsop50[ order( resSigGenesm9vsop50$log2FoldChange ), ] )
# strongest upregulation
tail( resSigGenesm9vsop50[ order( resSigGenesm9vsop50$log2FoldChange ), ] )








