#### Plot comparison of coverage across annotation types ####

# Load packages
library(ggplot2)

# Load coverage table
# coverage_table <- read.table(file = 'output/annotations/coverage/coverage_table.txt', header = TRUE)
coverage_table <- read.table(file = snakemake@input$coverage_table, header = TRUE)

# Subset
coverage_table <- subset(coverage_table, annotation_type == 'reference' |
                         annotation_type == 'liftover' |
                         annotation_type == 'combinedmakeronly' |
                         (annotation_type == 'stringtie' & diet == 'all') |
                         annotation_type == 'combinedall')

coverage_table$coverage_mb <- coverage_table$coverage / 1000000

# Reorder
coverage_table$annotation_type <- factor(coverage_table$annotation_type,
                                         levels = c('reference', 'liftover',
                                                    'combinedmakeronly', 'stringtie',
                                                    'combinedall'))

coverage_table$location <- factor(coverage_table$location,
                                  levels = c('vc2010', 'bristol', 'altadena'))

# Plot
# Labels
x_labels <- c('VC2010', 'N2', 'PS2025')

legend_labels <- c('VC2010 reference', 'Liftover',
                   'MAKER2',
                   'StringTie',
                   'Liftover plus new genes')

# Plot
tiff(snakemake@output$coverage_plot_tiff, width = 1600, height = 800, res = 300)
ggplot2::ggplot(data = coverage_table,
                mapping = aes(x = location, y = coverage_mb)) +
  geom_bar(stat = 'identity', aes(fill = annotation_type), position = position_dodge2(width = 0.9, preserve = "single")) +
  scale_fill_brewer(labels = legend_labels,  palette = 'Dark2', direction = -1) +
  theme(text = element_text(size = 10), axis.text.y=element_text(size = 6)) +
  xlab(element_blank()) +
  ylab('Coverage (Mb)') +
  theme_light() +
  theme(panel.grid = element_blank()) +
  theme(legend.key = element_rect(color = NA, fill = NA),
        legend.key.size = unit(0.5, "cm"),
        legend.position="bottom") +
  scale_y_continuous(expand = c(0, 0), labels = scales::comma, limits = c(0, 75)) +
  scale_x_discrete(labels = x_labels) +
  labs(fill = "Annotation type")
dev.off()

p <- last_plot()
# saveRDS(p, 'output/annotations/coverage/plots/coverage_plot_comparison.rds')
saveRDS(p, snakemake@output$coverage_plot_rds)
