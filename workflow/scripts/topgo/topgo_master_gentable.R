#### Combine gentables across different comparisons ####

# Load packages
library(tidyverse)

# load filepaths
# filepaths <- c('output/topgo/gentable/classicFisher/bristol/bristol_m9_op50.rds',
#                 'output/topgo/gentable/classicFisher/bristol/bristol_m9_pf.rds',
#                 'output/topgo/gentable/classicFisher/bristol/bristol_pf_op50.rds',
#                 'output/topgo/gentable/classicFisher/bristol/bristol_pf_as.rds',
#                 'output/topgo/gentable/classicFisher/altadena/altadena_m9_op50.rds',
#                 'output/topgo/gentable/classicFisher/altadena/altadena_m9_pf.rds',
#                 'output/topgo/gentable/classicFisher/altadena/altadena_pf_op50.rds',
#                 'output/topgo/gentable/classicFisher/altadena/altadena_pf_as.rds')
filepaths <- snakemake@input$topgo_gentable

# Combine the gentables for all the samples

master_gentable <- data.frame()

make_master_gentable <- function(filepath)
{
  # comparison
  name <- basename(filepath) %>%
    stringr::str_remove('.rds')

  gentable <- readRDS(filepath) %>%
    dplyr::mutate(name = name)

  master_gentable <<- rbind(master_gentable, gentable)
}

lapply(filepaths, make_master_gentable)

master_gentable <- cbind(master_gentable, str_split_fixed(master_gentable$name, '_', 2))
names(master_gentable)[names(master_gentable) == '1'] <- 'location'
names(master_gentable)[names(master_gentable) == '2'] <- 'comparison'

master_gentable <- mutate(master_gentable, classicFisher = stringr::str_replace(classicFisher, '< ', '')) %>%
  dplyr::mutate(classicFisher = as.numeric(classicFisher)) %>%
  dplyr::mutate(logpadj = -log(classicFisher))

saveRDS(master_gentable, file = snakemake@output$master_gentable)
