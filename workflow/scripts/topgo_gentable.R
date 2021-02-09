#### Generate table of most represented classes of differentially expressed genes ####

# Set working directory
setwd("~/OneDrive/Documents/Uni/III/Project/github/celegans-polyA/workflow")

# Load packages
library(topGO)

# Load GOdata
# load('output/topgo_data/bristol/bristol_as_op50.RData')
load(snakemake@input$topgo_data)

# Load topgo_results
# load('output/topgo_result/classic/fisher/bristol/bristol_as_op50.RData')
load(snakemake@input$topgo_result_classic_fisher)
topgo_result_classic_fisher <- topgo_result

# load('output/topgo_result/classic/ks/bristol/bristol_as_op50.RData')
load(snakemake@input$topgo_result_classic_ks)
topgo_result_classic_ks <- topgo_result

# load('output/topgo_result/elim/ks/bristol/bristol_as_op50.RData')
load(snakemake@input$topgo_result_elim_ks)
topgo_result_elim_ks <- topgo_result

all_res <- topGO::GenTable(GOdata,  classicFisher = topgo_result_classic_fisher,
                          classicKS = topgo_result_classic_ks,
                          elimKS = topgo_result_elim_ks,
                          orderBy = snakemake@wildcards$orderBy, 
                          ranksOf = "classicFisher", 
                          topNodes = 10)

save(all_res, file = snakemake@output$topgo_gentable)
