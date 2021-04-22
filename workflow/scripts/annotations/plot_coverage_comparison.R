#### Plot comparison of coverage across annotation types ####

# Set working directory
setwd("~/OneDrive/Documents/Uni/III/Project/github/celegans-polyA/workflow")

# Load packages
library(ggplot2)

# Load coverage table
# coverage_table <- read.table(file = 'output/annotations/coverage/coverage_table.txt', header = TRUE)
coverage_table <- read.table(file = snakemake@input$coverage_table, header = TRUE)

# Subset
coverage_table <- coverage_table[which(coverage_table$replicate != 'rep1_2_3' & 
                                         coverage_table$replicate != 'rep123' |
                                         coverage_table$annotation_type == 'liftover' &
                                         coverage_table$location == 'bristol'), ]

# Reorder
coverage_table$diet <- factor(coverage_table$diet, levels = rev(c("as", "bp", "hb101", "op50", "pf", "m9")))

# Labels
x_labels <- rev(c('Liftover', expression(italic('Acinetobacter\n  schindleri')),
              'Starvation'))

legend_labels <- c('Liftover', 'SOAPdenovo-trans', 'Stringtie', 'Trinity')

# Plot
tiff(snakemake@output$coverage_plot, width = 1600, height = 800, res = 300)
ggplot2::ggplot(data = coverage_table,
                mapping = aes(x = diet, y = coverage)) +
  geom_bar(stat = 'identity', aes(fill = annotation_type), position = position_dodge2(width = 0.9, preserve = "single")) +
  scale_fill_brewer(labels = legend_labels, palette = 'Set2', direction = -1) +
  coord_flip() +
  theme(text = element_text(size = 10), axis.text.y=element_text(size = 6)) +
  xlab(element_blank()) +
  ylab('Coverage (bp)') +
  theme_light() +
  scale_x_discrete(labels = x_labels) +
  scale_y_continuous(labels = scales::comma) +
  labs(fill = "Annotation type")
dev.off()

