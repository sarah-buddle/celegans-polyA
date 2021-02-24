#### Plot number of new genes by diet ####

# Set working directory
setwd("~/OneDrive/Documents/Uni/III/Project/github/celegans-polyA/workflow")

# Load packages
library(ggplot2)

# Load data
new_gene_counts <- read.table(file = 'output/annotations/new_genes/new_gene_counts.txt', 
                              header = TRUE)

# Labels
legend_labels <- c('Altadena', 'Bristol')

# Plot
ggplot2::ggplot(data = new_gene_counts, aes(x = diet, y = new_gene_count, fill = location)) +
  geom_bar(stat = 'identity', position = 'dodge2') +
  xlab(label = 'Diet') +
  ylab(label = 'Number of new genes') +
  theme_light() +
  theme(text = element_text(size = 14)) +
  labs(fill = 'Location') +
  scale_fill_brewer(labels = legend_labels)
  
