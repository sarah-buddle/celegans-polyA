#### Enrichment test in gene ontology analysis of differentially expressed genes ####

# Load packages
library(topGO)

# Load topgo_data
topgo_data <- readRDS(snakemake@input$topgo_data)

# Enrichment test
topgo_result <- topGO::runTest(topgo_data,
                                algorithm = snakemake@wildcards$algorithm,
                                statistic = snakemake@wildcards$statistic)
# Save result
saveRDS(topgo_result, file = snakemake@output$topgo_result)
