#### Enrichment test in gene ontology analysis of differentially expressed genes ####

# Set working directory
setwd("~/OneDrive/Documents/Uni/III/Project/github/celegans-polyA/workflow")

# Load packages
library(topGO)

# Load topgo_data
load(snakemake@input$topgo_data)

# Enrichment test
topgo_result <- topGO::runTest(GOdata,
                                algorithm = snakemake@wildcards$algorithm,
                                statistic = snakemake@wildcards$statistic)
# Save result
save(topgo_result, file = snakemake@output$topgo_result)