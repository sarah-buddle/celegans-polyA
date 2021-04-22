#### Gene ontology analysis of differentially expressed genes ####

# Load packages
library(topGO)

# Load mappings
geneID2GO <- topGO::readMappings(file = snakemake@input$wbid2go)

# List of all genes
geneUniverse <- names(geneID2GO)

# Load dds_res
dds_res <- readRDS(file = snakemake@input$dds_res)

# Create list of most differentially expressed genes

# tibble with the gene names and adjusted p values
p_values <- data.frame(rownames(dds_res), dds_res$padj)
colnames(p_values) <- c("gene_id", "padj")
p_values

# genes with padj < 0.05
topDiffGenes <- subset(p_values, padj < 0.05)
genesOfInterest <- as.character(topDiffGenes$gene_id)

# Tells topGO which genes in the gene universe are of interest
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

# Combine the above into topGOdata object
GOdata <- new("topGOdata",
              ontology="BP", # BP = biological process, MF = molecular function, CC = cellular component
              allGenes=geneList,
              annot = annFUN.gene2GO, # using user-defined mapping of genes to GO terms
              gene2GO = geneID2GO,
              nodeSize = 10) # Avoids problems with GP terms with few annotated genes being enriched due to statistical artefact

# Save topGO data object
saveRDS(GOdata, file = snakemake@output$topgo_data)
