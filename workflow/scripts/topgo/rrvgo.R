#### Order GO terms by similarity ####

# Load packages
library(rrvgo)

# Load data
master_gentable <- readRDS(snakemake@input$master_gentable)

# Order GO terms with by function

sim_matrix <- rrvgo::calculateSimMatrix(master_gentable$GO.ID,
                                       orgdb = 'org.Ce.eg.db',
                                       ont = 'BP',
                                       method = 'Rel')

reduced_terms <- reduceSimMatrix(sim_matrix,
                                threshold=0.7,
                                orgdb='org.Ce.eg.db')

colnames(reduced_terms)[which(names(reduced_terms) == 'go')] <- 'GO.ID'

saveRDS(reduced_terms, file = snakemake@output$rrvgo)
