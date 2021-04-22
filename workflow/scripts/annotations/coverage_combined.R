#### Calculare coverage of a combined annotation ####

# Load packages
library(GenomicRanges)
library(rtracklayer)

# Load annotation
annotation <- rtracklayer::import(snakemake@input$annotation)

# Reduce overlapping ranges
reduced_annotation <- GenomicRanges::reduce(annotation)

# Calculate coverage
coverage <- sum(width(reduced_annotation))

# Save output
saveRDS(coverage, file = snakemake@output$coverage)

# Write output to new row of coverage_table

# Forms sample name e.g. soap_bristol_as_rep1 or liftover_bristol
sample_name <- paste(snakemake@wildcards$source,
                       snakemake@wildcards$location, 'all_rep123', sep = "_")

new_row <- c(sample_name, snakemake@wildcards$source,
               snakemake@wildcards$location, 'all', 'rep123', coverage)


# Making or updating coverage_table.txt
# Checks whether coverage_table.txt file already exists
if(file.exists('output/annotations/coverage/coverage_table.txt') == TRUE) {

  # import existing table
  coverage_table <- read.table('output/annotations/coverage/coverage_table.txt',
                               header = TRUE)

  # removes row corresponding to current sample if necessary
  coverage_table <- coverage_table[coverage_table$sample_name != sample_name, ]

  # adds new row
  coverage_table <- rbind(coverage_table, new_row)

  # sort alphabetically by sample name
  coverage_table <- coverage_table[order(coverage_table$sample_name), ]

  # write  to coverage_table.txt
  write.table(coverage_table, file = 'output/annotations/coverage/coverage_table.txt',
              row.names = FALSE)

  # if no coverage_table.txt file
}else{
  # column headings
  headers <- c('sample_name', 'annotation_type', 'location', 'diet', 'replicate', 'coverage')

  # Makes data frame from headings nd new row
  coverage_table <- t(data.frame(headers, new_row))

  # Write to coverage_table.txt
  write.table(coverage_table, file = 'output/annotations/coverage/coverage_table.txt',
              col.names = FALSE, row.names = FALSE)
}
