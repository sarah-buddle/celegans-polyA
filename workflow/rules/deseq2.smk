LOCATIONS = ['bristol','altadena']
DIETS = ['as','bp','hb101','m9','op50','pf']
REPLICATES = ['rep1','rep2','rep3']

rule deseq2:
    ''' Creates DESeq2 objects for use in further analysis '''
    input:
        script='scripts/deseq2/deseq2.R',
        samples='scripts/deseq2/htseqcount_samples_full.csv',
        htseq_count=expand('from_cluster/htseq_count_export/{{location}}_{diet}_{replicate}_counts.txt', \
        diet=DIETS, replicate=REPLICATES)
    output:
        rlog_dds='output/deseq2_data/{location}/{location}_rlog_dds.RData',
        full_dds='output/deseq2_data/{location}/{location}_full_dds.RData'
    conda:
        '../envs/conda/bioconductor-deseq2=1.30.0.yaml'
    script:
        '../scripts/deseq2/deseq2.R'

rule deseq2_PCA:
    ''' PCA Plot to summarise differential expression analysis '''
    input:
        script='scripts/deseq2/deseq2_pca.R',
        rlog_dds='output/deseq2_data/{location}/{location}_rlog_dds.RData'
    output:
        pca_plot='output/deseq2_plots/{location}/{location}_deseq2_pca.tiff'
    conda:
        '../envs/conda/bioconductor-deseq2=1.30.0_r-ggplot2=3.3.1.yaml'
    script:
        '../scripts/deseq2/deseq2_pca.R'

'''
snakemake --cores 1 --use-conda output/deseq2_plots/altadena/altadena_deseq2_pca.tiff
'''

rule deseq2_no_m9:
    ''' Creates DESeq2 object without m9/starvation diet treatment '''
    input:
        script='scripts/deseq2/deseq2_no_m9.R',
        samples='scripts/htseqcount_samples_full.csv',
        htseq_count=expand('from_cluster/htseq_count_export/{{location}}_{diet}_{replicate}_counts.txt', \
        diet=DIETS, replicate=REPLICATES)
    output:
        rlog_dds='output/deseq2_data/{location}/{location}_no_m9_rlog_dds.RData'
    conda:
        '../envs/conda/bioconductor-deseq2=1.30.0.yaml'
    script:
        '../scripts/deseq2/deseq2_no_m9.R'

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
    ''' Plot percentage genes expressed for single location '''
    input:
        script='scripts/deseq2/plot_percentage_genes_expressed.R',
        full_dds='output/deseq2_data/{location}/{location}_full_dds.RData',
        samples='scripts/htseqcount_samples_full.csv'
    output:
        percentage_plot='output/deseq2_plots/{location}/{location}_percentage_plot.png'
    conda:
        '../envs/conda/bioconductor-deseq2=1.30.0_r-ggplot2=3.3.1.yaml'
    script:
        '../scripts/deseq2/plot_percentage_genes_expressed.R'

'''
snakemake --cores 1 --use-conda \
output/deseq2_plots/altadena/altadena_percentage_plot.png
'''

rule count_percentage_genes_expressed:
    ''' Count percentage of total genes expressed '''
    input:
        script='scripts/deseq2/count_percentage_genes_expressed.R',
        full_dds='output/deseq2_data/{location}/{location}_full_dds.RData',
        samples='scripts/htseqcount_samples_full.csv'
    output:
        expressed_genes='output/expressed_genes/{location}/{location}_expressed_genes.RData'
    conda:
        '../envs/conda/bioconductor-deseq2=1.30.0.yaml'
    script:
        '../scripts/deseq2/count_percentage_genes_expressed.R'

'''
snakemake --cores 1 --use-conda \
output/expressed_genes/bristol/bristol_expressed_genes.RData
'''

rule plot_percentage_expression_all:
    input:
        script='scripts/plot_percentage_expression_all.R',
        full_dds_altadena='output/expressed_genes/altadena/altadena_expressed_genes.RData',
        full_dds_bristol='output/expressed_genes/bristol/bristol_expressed_genes.RData'
    output:
        percentage_plot='output/deseq2_plots/all/percentage_plot_all.tiff'
    conda:
        '../envs/conda/bioconductor-deseq2=1.30.0_r-ggplot2=3.3.1.yaml'
    script:
        '../scripts/plot_percentage_expression_all.R'

'''
snakemake --cores 1 --use-conda \
output/deseq2_plots/all/percentage_plot_all.tiff
'''

rule de_analysis:
    input:
        script='scripts/deseq2/de_analysis.R',
        full_dds='output/deseq2_data/{location}/{location}_full_dds.RData'
    output:
        dds_res='output/deseq2_results/{location}/{location}_{diet1}_{diet2}.RData'
    conda:
        '../envs/conda/bioconductor-deseq2=1.30.0.yaml'
    script:
        '../scripts/deseq2/de_analysis.R'

'''
snakemake --cores 1 --use-conda output/deseq2_data/bristol/bristol_as_op50.RData
'''

rule deseq2_all_count_matrix:
    input:
        script='scripts/deseq2/deseq2_all_count_matrix.R',
        htseq_count=expand('from_cluster/htseq_count_export/{location}_{diet}_{replicate}_counts.txt', \
        location = LOCATIONS, diet=DIETS, replicate=REPLICATES)
    output:
        count_matrix='output/deseq2_data/all/count_matrix.rds'
    conda:
        '../envs/conda/r-tidyverse=1.2.1.yaml'
    script:
        '../scripts/deseq2/deseq2_all_count_matrix.R'

'''
snakemake --cores 1 --use-conda -R \
output/deseq2_data/all/count_matrix.rds
'''

rule deseq2_all:
    input:
        script='scripts/deseq2/deseq2_all.R',
        samples='scripts/deseq2/htseqcount_samples_full.csv',
        count_matrix='output/deseq2_data/all/count_matrix.rds'
    output:
        vst_dds='output/deseq2_data/all/all_vst_dds.rds',
        full_dds='output/deseq2_data/all/all_full_dds.rds'
    conda:
        '../envs/conda/bioconductor-deseq2=1.30.0.yaml'
    script:
        '../scripts/deseq2/deseq2_all.R'

'''
snakemake --cores 1 --use-conda -R \
output/deseq2_data/all/all_vst_dds.rds
'''

rule deseq2_pca_all:
    ''' PCA Plot to summarise differential expression analysis '''
    input:
        script='scripts/deseq2/all_deseq2_pca.R',
        vst_dds='output/deseq2_data/all/all_vst_dds.rds'
    output:
        pca_plot='output/deseq2_plots/all/all_deseq2_pca.tiff'
    conda:
        '../envs/conda/bioconductor-deseq2=1.30.0_r-ggplot2=3.3.1.yaml'
    script:
        '../scripts/deseq2/all_deseq2_pca.R'

'''
snakemake --cores 1 --use-conda -R \
output/deseq2_plots/all/all_deseq2_pca.tiff
'''
