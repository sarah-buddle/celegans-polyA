LOCATIONS = ['bristol','altadena']
DIETS = ['as','bp','hb101','m9','op50','pf']
REPLICATES = ['rep1','rep2','rep3']

# rule all:
#     input:
#         'output/deseq2/{location}_percentage_plot.png'

rule deseq2:
    ''' Creates DESeq2 objects for use in further analysis '''
    input:
        script='scripts/deseq2.R',
        samples='scripts/htseqcount_samples_full.csv',
        htseq_count=expand('from_cluster/htseq_count_export/{{location}}_{diet}_{replicate}_counts.txt', \
        diet=DIETS, replicate=REPLICATES)
    output:
        rlog_dds='output/deseq2_data/{location}/{location}_rlog_dds.RData',
        full_dds='output/deseq2_data/{location}/{location}_full_dds.RData'
    conda:
        '../envs/conda/bioconductor-deseq2=1.30.0.yaml'
    script:
        '../scripts/deseq2.R'

'''
snakemake --cores 1 --use-conda output/deseq2_data/bristol/bristol_full_dds.RData
'''

rule deseq2_PCA:
    ''' PCA Plot to summarise differential expression analysis '''
    input:
        script='scripts/deseq2_pca.R',
        rlog_dds='output/deseq2_data/{location}/{location}_rlog_dds.RData'
    output:
        pca_plot='output/deseq2_plots/{location}/{location}_deseq2_pca.png'
    conda:
        '../envs/conda/bioconductor-deseq2=1.30.0_r-ggplot2=3.3.1.yaml'
    script:
        '../scripts/deseq2_pca.R'

'''
snakemake --cores 1 --use-conda output/deseq2_plots/bristol/bristol_deseq2_pca.png
'''

rule deseq2_no_m9:
    ''' Creates DESeq2 object without m9/starvation diet treatment '''
    input:
        script='scripts/deseq2.R',
        samples='scripts/htseqcount_samples_full.csv',
        htseq_count=expand('from_cluster/htseq_count_export/{{location}}_{diet}_{replicate}_counts.txt', \
        diet=DIETS, replicate=REPLICATES)
    output:
        rlog_dds='output/deseq2_data/{location}/{location}_no_m9_rlog_dds.RData'
    conda:
        '../envs/conda/bioconductor-deseq2=1.30.0.yaml'
    script:
        '../scripts/deseq2_no_m9.R'

'''
snakemake --cores 1 --use-conda output/deseq2/bristol_no_m9_rlog_dds.RData
'''

rule deseq2_PCA_no_m9:
    input:
        script='scripts/deseq2_pca.R',
        rlog_dds='output/deseq2_data/{location}/{location}_no_m9_rlog_dds.RData'
    output:
        pca_plot='output/deseq2_plots/{location}/{location}_deseq2_pca_no_m9.png'
    conda:
        '../envs/conda/bioconductor-deseq2=1.30.0_r-ggplot2=3.3.1.yaml'
    script:
        '../scripts/deseq2_pca.R'

'''
snakemake --cores 1 --use-conda output/deseq2/bristol_deseq2_pca_no_m9.png
'''

# rule count_genes:
#     ''' Counts unique WB IDs in liftover annotation '''
#     input:
#         script='scripts/count_genes.R',
#         liftover='from_cluster/liftover_annotations/{location}/liftover_{location}.gtf'
#     output:
#         gene_count='output/gene_counts/{location}/{location}_gene_count.RData'
#     conda:
#         '../envs/conda/bioconductor-rtracklayer=1.50.0.yaml'
#     script:
#         '../scripts/count_genes.R'
#
# '''
# snakemake --cores 1 --use-conda output/gene_counts/bristol_gene_count.RData
# '''

rule deseq2_percentage_genes_expressed:
    input:
        script='scripts/deseq2_percentage_genes_expressed.R',
        full_dds='output/deseq2_data/{location}/{location}_full_dds.RData',
        samples='scripts/htseqcount_samples_full.csv'
    output:
        percentage_plot='output/deseq2_plots/{location}/{location}_percentage_plot.png'
    conda:
        '../envs/conda/bioconductor-deseq2=1.30.0_r-ggplot2=3.3.1.yaml'
    script:
        '../scripts/deseq2_percentage_genes_expressed.R'

'''
snakemake --cores 1 --use-conda \
output/deseq2_plots/altadena/altadena_percentage_plot.png
'''

rule de_analysis:
    input:
        script='scripts/de_analysis.R',
        full_dds='output/deseq2_data/{location}/{location}_full_dds.RData'
    output:
        dds_res='output/deseq2_data/{location}/{location}_{diet1}_{diet2}.RData'
    conda:
        '../envs/conda/bioconductor-deseq2=1.30.0.yaml'
    script:
        '../scripts/de_analysis.R'

'''
snakemake --cores 1 --use-conda output/deseq2_data/bristol/bristol_as_op50.RData
'''
