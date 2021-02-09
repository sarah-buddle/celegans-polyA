#### Plot Coverage ####

# Set working directory
setwd("~/OneDrive/Documents/Uni/III/Project/github/celegans-polyA/workflow")

# Load packages
library(ggplot2)

# Load coverage table
# coverage_table <- read.table(file = 'output/annotations/coverage/coverage_table.txt', header = TRUE)
coverage_table <- read.table(file = snakemake@input$coverage_table, header = TRUE)

# Subset
coverage_table <- coverage_table[which(coverage_table$annotation_type == 'soap',
                                       coverage_table$location == 'bristol', 
                                       coverage_table$replicate == 'rep1_2_3'), ]

# Reorder
coverage_table$diet <- factor(coverage_table$diet, levels = rev(c("as", "bp", "hb101", "op50", "pf", "m9")))

# Labels
x_labels <- c(expression(italic("Acinetobacter schindleri")),
              expression(italic("Bacillus pumilus")),
              expression(paste(italic("Escherichia coli"), " strain HB101")),
              expression(paste(italic("Escherichia coli"), " strain OP50")),
              expression(italic("Pseudomonas fragi")),
              "Starvation")

# Plot
png(snakemake@output$coverage_plot)
ggplot2::ggplot(data = coverage_table,
                mapping = aes(x = diet, y = coverage)) +
  geom_bar(stat = 'identity', aes(fill = diet)) +
  scale_fill_brewer(palette = 'Set2', direction = -1) +
  coord_flip() +
  scale_x_discrete(labels = rev(x_labels)) +
  theme(legend.position = "none") +
  theme(text = element_text(size=14)) +
  xlab('Diet') +
  ylab('coverage')
dev.off()

