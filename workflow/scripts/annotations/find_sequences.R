# Set working directory
setwd("~/OneDrive/Documents/Uni/III/Project/github/celegans-polyA/workflow")

# Load packages
library(seqinr)

# Load file paths for genomes
# genome_paths <- c('input/genomes/altadena_genome_chr/altadena_genome_chrI.fa', 
#                   'input/genomes/altadena_genome_chr/altadena_genome_chrII.fa', 
#                   'input/genomes/altadena_genome_chr/altadena_genome_chrIII.fa', 
#                   'input/genomes/altadena_genome_chr/altadena_genome_chrIV.fa', 
#                   'input/genomes/altadena_genome_chr/altadena_genome_chrV.fa', 
#                   'input/genomes/altadena_genome_chr/altadena_genome_chrX.fa')
genome_paths <- snakemake@input$genome

# Load genome
genome <- lapply(genome_paths, seqinr::read.fasta, as.string = TRUE)
names(genome) <- c('I','II','III','IV','V','X')

# Load granges_new_genes
# new_genes <- data.frame(readRDS('output/annotations/granges_new_genes_merged/altadena/altadena.rds'))
new_genes <- data.frame(readRDS(snakemake@input$granges_new_genes))

# Extract sequence using ranges
# extract_sequence <- function(rownumber)
# {
#   start <- new_genes$start[rownumber]
#   end <- new_genes$end[rownumber]
#   chromosome <- new_genes$seqnames[rownumber]
#   sequence <- substr(genome[chromosome], start, end)
#   seqinr::write.fasta(sequence, as.string = TRUE, 
#                       name = paste(chromosome, '_', rownumber, sep = ''), 
#                       file.out = paste('output/annotations/new_gene_sequences/',
#                                        snakemake@wildcards$location,
#                                        chromosome, '_', rownumber, '.fa', 
#                                        sep = ''))
# }

extract_sequence <- function(rownumber)
{
  start <- new_genes$start[rownumber]
  end <- new_genes$end[rownumber]
  chromosome <- new_genes$seqnames[rownumber]
  sequence <- substr(genome[chromosome], start, end)
  seqinr::write.fasta(sequence, as.string = TRUE, 
                      name = paste(chromosome, '_', rownumber, sep = ''), 
                      file.out = snakemake@output$sequences, open = "a")
}

lapply(1:nrow(new_genes), extract_sequence)
