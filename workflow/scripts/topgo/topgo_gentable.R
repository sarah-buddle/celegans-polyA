#### Generate table of most represented classes of differentially expressed genes ####

# Load packages
library(topGO)

# Load GOdata
#load('output/topgo_data/altadena/altadena_as_op50.RData')
topgo_data <- readRDS(snakemake@input$topgo_data)

# Load topgo_results
#load('output/topgo_result/classic/fisher/altadena/altadena_as_op50.RData')
topgo_result_classic_fisher <- readRDS(snakemake@input$topgo_result_classic_fisher)

#load('output/topgo_result/classic/ks/altadena/altadena_as_op50.RData')
topgo_result_classic_ks <- readRDS(snakemake@input$topgo_result_classic_ks)

#load('output/topgo_result/elim/ks/altadena/altadena_as_op50.RData')
topgo_result_elim_ks <-readRDS(snakemake@input$topgo_result_elim_ks)

all_res <- topGO::GenTable(topgo_data,  classicFisher = topgo_result_classic_fisher,
                          classicKS = topgo_result_classic_ks,
                          elimKS = topgo_result_elim_ks,
                          #orderBy = 'classicFisher',
                          orderBy = snakemake@wildcards$orderBy,
                          ranksOf = "classicFisher",
                          topNodes = 1000)

saveRDS(all_res, file = snakemake@output$topgo_gentable)
