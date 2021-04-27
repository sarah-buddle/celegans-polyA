# Finds sequence of new genes

# Load packages
library(seqinr)

# Load file paths for genomes
# genome_paths <- c('input/genomes/bristol_genome_chr/bristol_genome_chrI.fa',
#                   'input/genomes/bristol_genome_chr/bristol_genome_chrII.fa',
#                   'input/genomes/bristol_genome_chr/bristol_genome_chrIII.fa',
#                   'input/genomes/bristol_genome_chr/bristol_genome_chrIV.fa',
#                   'input/genomes/bristol_genome_chr/bristol_genome_chrV.fa',
#                   'input/genomes/bristol_genome_chr/bristol_genome_chrX.fa')
genome_paths <- snakemake@input$genome

# Load genome
genome <- lapply(genome_paths, seqinr::read.fasta, as.string = TRUE)
names(genome) <- c('I','II','III','IV','V','X')

# Load granges_new_genes
# new_genes <- data.frame(readRDS('output/annotations/combined_annotations/combinedall/bristol/allnewgenes_bristol.rds'))
new_genes <- data.frame(readRDS(snakemake@input$allnewgenes))

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
  gene_id <- new_genes$gene_id[rownumber]
  start <- new_genes$start[rownumber]
  end <- new_genes$end[rownumber]
  chromosome <- new_genes$seqnames[rownumber]
  sequence <- substr(genome[chromosome], start, end)
  seqinr::write.fasta(sequence, as.string = TRUE,
                      # name = paste('bristol', gene_id, sep = '_'),
                      name = paste(snakemake@wildcards$location, gene_id, sep = '_'),
                      nbchar = 50000,
                      # file.out = 'output/blast/new_gene_sequences/bristol/bristol.fa', open = "a")
                      file.out = snakemake@output$sequences, open = "a")
}

lapply(1:nrow(new_genes), extract_sequence)
