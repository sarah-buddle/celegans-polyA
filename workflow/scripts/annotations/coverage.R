#### Calculate coverage of the annotations #### 

# Set working directory
setwd("~/OneDrive/Documents/Uni/III/Project/github/celegans-polyA/workflow")

# Load packages
library(GenomicRanges)

# Load annotation
load(file = snakemake@input$granges)

# Reduce overlapping ranges
reduced_annotation <- GenomicRanges::reduce(annotation)

# Calculate coverage
coverage <- sum(width(reduced_annotation))

# Save output
save(coverage, file = snakemake@output$coverage)

# Write output to new row of coverage_table

# Forms sample name e.g. soap_bristol_as_rep1 or liftover_bristol
if(snakemake@wildcards$annotation_type == 'liftover'){
  sample_name <- paste('liftover_', snakemake@wildcards$location, sep = "")
  
  # New row
  # new_row <- c(sample_name, 'soap', 'bristol', 'as', 'rep2', 100)
  new_row <- c(sample_name, 'liftover', snakemake@wildcards$location, NA, NA, coverage)
  
}else{
  #sample_name <- 'soap_bristol_as_rep2'
  sample_name <- paste(snakemake@wildcards$annotation_type, '_',
                        snakemake@wildcards$location, '_', 
                        snakemake@wildcards$diet, '_', 
                        snakemake@wildcards$replicate, sep = "")
  
  # New row
  # new_row <- c(sample_name, 'soap', 'bristol', 'as', 'rep2', 100)
  new_row <- c(sample_name, snakemake@wildcards$annotation_type, 
               snakemake@wildcards$location, snakemake@wildcards$diet, 
               snakemake@wildcards$replicate, coverage)
}

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
coverage_table <- read.table('output/annotations/coverage/coverage_table.txt', header = TRUE)
coverage_table <- coverage_table[order(coverage_table$sample_name), ]


