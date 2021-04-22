# Find sequences of new genes from their names

# Set working directory
setwd('~/OneDrive/Documents/Uni/III/Project/github/celegans-polyA/workflow')

# Load packages
library(seqinr)

# Load file
no_hits_names <- read.table('output/annotations/blast_no_hits/altadena/altadena_blast_no_hits_names.txt')
no_hits_names <- read.table(snakemake@input$no_hits_names)
no_hits_names <- as.list(no_hits_names$V1)

# Load sequences
sequences <- seqinr::read.fasta('output/annotations/new_gene_sequences/altadena/altadena.fa')
sequences <- seqinr::read.fasta(snakemake@input$sequences)

no_hits <- sequences[c(which(names(sequences) %in% no_hits_names))]

seqinr::write.fasta(no_hits, file = 'output/annotations/blast_no_hits/altadena/altadena_blast_no_hits.fa',
                    names = names(no_hits), nbchar = 50000)
seqinr::write.fasta(no_hits, file = snakemake@output$no_hits,
                    names = names(no_hits), nbchar = 50000)
