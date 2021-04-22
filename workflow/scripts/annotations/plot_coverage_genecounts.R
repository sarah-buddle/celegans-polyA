#### Combine coverage and gene counts plots ####

# Load packages
library(ggplot2)
library(gridExtra)
library(ggpubr)

# Load plots
# coverage <- readRDS('output/annotations/coverage/plots/coverage_plot_comparison.rds')
coverage <- readRDS(snakemake@input$coverage_plot)

#genes <- readRDS('output/annotations/gene_counts/gene_counts_plot.rds')
genes <- readRDS(snakemake@input$gene_plot)

# Extract legend
legend <- ggpubr::get_legend(coverage)

# Remove legends
coverage <- coverage + theme(legend.position = 'none')
genes <- genes + theme(legend.position = 'none')

# Plot
tiff(snakemake@output$coverage_genecounts_plot, width = 2000, height = 1000, res = 300)
gridExtra::grid.arrange(coverage, genes, legend, ncol=2, nrow = 2,
                        layout_matrix = rbind(c(1,2), c(3,3)),
                        widths = c(2.7, 2.7), heights = c(2.5, 0.2))
dev.off()
