#### Plot expression and annotation type of new genes ####

# Load packages
library(ggplot2)

# Load data
# counts_bristol <- readRDS('output/annotations/new_gene_expression/data/bristol/new_gene_expression_bristol.rds')
counts_bristol <- readRDS(snakemake@input$expression_counts_bristol)
counts_bristol$location <-  'bristol'

# counts_altadena <- readRDS('output/annotations/new_gene_expression/data/altadena/new_gene_expression_altadena.rds')
counts_altadena <- readRDS(snakemake@input$expression_counts_altadena)
counts_altadena$location <- 'altadena'

data <- rbind(counts_bristol, counts_altadena)

# Reorder
data$annotation_type <- factor(data$annotation_type,
                               levels = c('maker', 'StringTie', 'maker_stringtie'))

data$genes <- factor(data$genes, levels = c('not_expressed', 'expressed'))

data$location <- factor(data$location, levels = c('bristol', 'altadena'))

# Labels

x_labels <- c('MAKER2\nonly', 'StringTie\nonly', 'MAKER2 and\nStringTie')
legend_labels <- c('Very low expression', 'Higher expression')

locations <- c('N2', 'PS2025')
names(locations) <- c('bristol', 'altadena')

# tiff('output/annotations/new_gene_expression/plots/new_gene_expression_plot.tiff', width = 1600, height = 800, res = 300)
tiff(snakemake@output$new_gene_expression_plot, width = 1600, height = 800, res = 300)
ggplot2::ggplot(data, mapping = aes(x = annotation_type, y = number_of_genes, fill = genes)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~ location, labeller = labeller(location = locations)) +
  scale_fill_brewer(labels = legend_labels) +
  xlab('Annotation type') +
  ylab('Number of predicted new genes') +
  theme_light() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key = element_rect(color = NA, fill = NA),
        legend.key.size = unit(0.5, "cm"),
        legend.position = c(0.14, 0.84)) +
  scale_x_discrete(labels = x_labels) +
  scale_y_continuous(expand = c(0, 0), labels = scales::comma, limits = c(0, 400))
dev.off()
