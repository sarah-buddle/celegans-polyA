# Produce phylogenetic tree of species used in blast analysis

# Set working directory
setwd('~/OneDrive/Documents/Uni/III/Project/github/celegans-polyA/workflow')

# Load packages
library(ggtree)

# Load data (from Cristian)
# tree <- ggtree::read.tree('input/tree/SpeciesTree_rooted_node_labels.txt')
tree <- ggtree::read.tree(snakemake@input$tree)

genus <- c('Trichuris ', 'Onchocerca ', 'Brugia ', 'Panagrellus ', 'Strongyloides ', 'Pristionchus ', 'Oscheius ', 'Caenorhabditis ', 'Caenorhabditis ', 'Caenorhabditis ', 'Caenorhabditis ', 'Caenorhabditis ', 'Caenorhabditis ', 'Caenorhabditis ', 'Caenorhabditis ', 'Caenorhabditis ', 'Caenorhabditis ', 'Caenorhabditis ')
species <- c('muris', 'volvulus', 'malayi', 'redivivus', 'ratti', 'pacificus', 'tipulae', 'angaria', 'japonica', 'inopinata', 'elegans', 'tropicalis', 'brenneri', 'sinica', 'briggsae', 'nigoni', 'latens', 'remanei')
clade <- c('Trichinellida', 'Spirurina', 'Spirurina', 'Tylenchina', 'Tylenchina', 'Rhabditina', 'Rhabditina', 'Rhabditina', 'Rhabditina', 'Rhabditina', 'Rhabditina', 'Rhabditina', 'Rhabditina', 'Rhabditina', 'Rhabditina', 'Rhabditina', 'Rhabditina', 'Rhabditina')
df <- data.frame(label = tree$tip.label, genus = genus,
                 species = species, clade = clade)

# tiff('output/blast/venn/bristol/venn_bristol.tiff', width = 1600, height = 800, res = 300)
tiff(snakemake@output$plot, width = 900, height = 1600, res = 300)
p <- ggtree::ggtree(tree) %<+% df +
  ggtree::xlim(0, 1.5) +
  ggtree::geom_tiplab(ggplot2::aes(label = paste0('italic(', genus, ')~italic(', species, ')')), parse = TRUE)
ggtree::viewClade(p, MRCA(p, 'celegans', 'cnigoni'))
dev.off()
