#### Comparing numbers of expressed genes between diet treatments ####

# Load packages
library(DESeq2)
library(ggplot2)

# Import sample list
samples <- as.data.frame(read.csv(snakemake@input$samples))

# Set location
samples <- subset(samples, location == snakemake@wildcards$location)

# Load appropriate DESeq object made by deseq scrip
full_dds <- readRDS(snakemake@input$full_dds)

# convert DESeq object into data frame
full_dds_df <- as.data.frame(assay(full_dds))
full_dds_df$gene_id <- rownames(full_dds_df)
rownames(full_dds_df) <- NULL

# Count total unique gene IDs
total_genes <- length(unique(full_dds_df$gene_id))

# function to count genes with reads above a certain threshold
count_expressed_genes <- function(df, sample, threshold)
{
  output <- sum(df[, sample] > threshold)
  return(output)
}

# apply the function to each sample in our data frame
sample_names <- samples$sample_name
expressed_gene_count <- lapply(sample_names, count_expressed_genes,
                               df = full_dds_df, threshold = 1)

# make data frame of results
expressed_genes <- do.call(rbind, Map(data.frame, sample_name = sample_names,
                                      expressed_gene_count=expressed_gene_count,
                                      location=samples$location,
                                      diet= samples$diet))
rownames(expressed_genes) <- NULL

# Calculate percentage of total genes expressed
expressed_genes$percentage_expressed_genes <-
  expressed_genes$expressed_gene_count*100/total_genes

# Plot percentages

# Legend
# Reorder
expressed_genes$diet <- factor(expressed_genes$diet,
                               levels = c("as", "bp", "hb101", "op50", "pf", "m9"))
# labels
legend_labels <- c(expression(italic("Acinetobacter schindleri")),
                    expression(italic("Bacillus pumilus")),
                    expression(paste(italic("Escherichia coli"), " strain HB101")),
                    expression(paste(italic("Escherichia coli"), " strain OP50")),
                    expression(italic("Pseudomonas fragi")),
                    "Starvation")

# x axis

# Reorder
location <- snakemake@wildcards$location

levels <- rev(c(paste(location, "_as_rep1", sep = ""), paste(location, "_as_rep2", sep = ""),
                paste(location, "_as_rep3", sep = ""), paste(location, "_bp_rep1", sep = ""),
                paste(location, "_bp_rep2", sep = ""), paste(location, "_bp_rep3", sep = ""),
                paste(location, "_hb101_rep1", sep = ""), paste(location, "_hb101_rep2", sep = ""),
                paste(location, "_hb101_rep3", sep = ""), paste(location, "_op50_rep1", sep = ""),
                paste(location, "_op50_rep2", sep = ""), paste(location, "_op50_rep3", sep = ""),
                paste(location, "_pf_rep1", sep = ""), paste(location, "_pf_rep2", sep = ""),
                paste(location, "_pf_rep3", sep = ""), paste(location, "_m9_rep1", sep = ""),
                paste(location, "_m9_rep2", sep = ""), paste(location, "_m9_rep3", sep = "")))

# levels <- rev(c("bristol_as_rep1", "bristol_as_rep2", "bristol_as_rep3",
#                 "bristol_bp_rep1", "bristol_bp_rep2", "bristol_bp_rep3",
#                 "bristol_hb101_rep1", "bristol_hb101_rep2", "bristol_hb101_rep3",
#                 "bristol_op50_rep1", "bristol_op50_rep2", "bristol_op50_rep3",
#                 "bristol_pf_rep1", "bristol_pf_rep2", "bristol_pf_rep3",
#                 "bristol_m9_rep1", "bristol_m9_rep2", "bristol_m9_rep3"))

expressed_genes$sample_name <- factor(expressed_genes$sample_name,
                                      levels = levels)

# Labels
x_axis_labels = rev(c("As rep1", "As rep2", "As rep3", "Bp rep1", "Bp rep2", "Bp rep3",
                "E.coli HB101 rep1", "E.coli HB101 rep2", "E.coli HB101 rep3",
                "E.coli OP50 rep1", "E.coli OP50 rep2", "E.coli OP50 rep3",
                "Pf rep1", "Pf rep2", "Pf rep3",
                "Starvation rep1", "Starvation rep2", "Starvation rep3"))

# Plot
tiff(snakemake@output$percentage_plot)
ggplot(data = expressed_genes, aes(x = sample_name,
                                   y = percentage_expressed_genes, fill = diet)) +
  geom_bar(stat = "identity") +
  theme(axis.title.y = element_blank()) +
  ylab("Percentage of total genes expressed") +
  coord_flip() +
  theme(text = element_text(size=14)) +
  labs(fill = "Diet") +
  scale_x_discrete(labels = x_axis_labels) +
  scale_fill_brewer(labels = legend_labels,
                    palette = "Set2")
dev.off()
