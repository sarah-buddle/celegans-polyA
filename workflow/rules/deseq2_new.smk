LOCATIONS = ['bristol','altadena']
DIETS = ['as','bp','hb101','m9','op50','pf']
REPLICATES = ['rep1','rep2','rep3']

rule deseq2_new:
    ''' Creates DESeq2 objects for use in further analysis '''
    input:
        script='scripts/deseq2/deseq2.R',
        samples='scripts/deseq2/htseqcount_samples_full.csv',
        htseq_count=expand('from_cluster/htseq_count_export_new/{{location}}_{diet}_{replicate}_counts.txt', \
        diet=DIETS, replicate=REPLICATES)
    output:
        rlog_dds='output/deseq2_new/data/{location}/{location}_rlog_dds.rds',
        full_dds='output/deseq2_new/data/{location}/{location}_full_dds.rds'
    conda:
        '../envs/conda/bioconductor-deseq2=1.30.0.yaml'
    script:
        '../scripts/deseq2/deseq2.R'

rule deseq2_PCA_new:
    ''' PCA Plot to summarise differential expression analysis '''
    input:
        script='scripts/deseq2/deseq2_pca.R',
        rlog_dds='output/deseq2_new/data/{location}/{location}_rlog_dds.rds'
    output:
        pca_plot='output/deseq2_new/plots/{location}/{location}_deseq2_pca.tiff'
    conda:
        '../envs/conda/bioconductor-deseq2=1.30.0_r-ggplot2=3.3.1.yaml'
    script:
        '../scripts/deseq2/deseq2_pca.R'

'''
snakemake --cores 1 --use-conda \
output/deseq2_new/plots/altadena/altadena_deseq2_pca.tiff \
output/deseq2_new/plots/bristol/bristol_deseq2_pca.tiff
'''

rule deseq2_no_m9_new:
    ''' Creates DESeq2 object without m9/starvation diet treatment '''
    input:
        script='scripts/deseq2/deseq2_no_m9.R',
        samples='scripts/htseqcount_samples_full.csv',
        htseq_count=expand('from_cluster/htseq_count_export_new/{{location}}_{diet}_{replicate}_counts.txt', \
        diet=DIETS, replicate=REPLICATES)
    output:
        rlog_dds='output/deseq2_new/data/{location}/{location}_no_m9_rlog_dds.rdata'
    conda:
        '../envs/conda/bioconductor-deseq2=1.30.0.yaml'
    script:
        '../scripts/deseq2/deseq2_no_m9.R'


rule deseq2_PCA_no_m9_new:
    input:
        script='scripts/deseq2_pca.R',
        rlog_dds='output/deseq2_new/data/{location}/{location}_no_m9_rlog_dds.rds'
    output:
        pca_plot='output/deseq2_plots/{location}/{location}_deseq2_pca_no_m9.tiff'
    conda:
        '../envs/conda/bioconductor-deseq2=1.30.0_r-ggplot2=3.3.1.yaml'
    script:
        '../scripts/deseq2_pca.R'


rule deseq2_percentage_genes_expressed_new:
    ''' Plot percentage genes expressed for single location '''
    input:
        script='scripts/deseq2/plot_percentage_genes_expressed.R',
        full_dds='output/deseq2_new/data/{location}/{location}_full_dds.rds',
        samples='scripts/htseqcount_samples_full.csv'
    output:
        percentage_plot='output/deseq2_new/plots/{location}/{location}_percentage_plot.tiff'
    conda:
        '../envs/conda/bioconductor-deseq2=1.30.0_r-ggplot2=3.3.1.yaml'
    script:
        '../scripts/deseq2/plot_percentage_genes_expressed.R'

'''
snakemake --cores 1 --use-conda \
output/deseq2_plots/altadena/altadena_percentage_plot.png
'''

rule count_percentage_genes_expressed_new:
    ''' Count percentage of total genes expressed '''
    input:
        script='scripts/deseq2/count_percentage_genes_expressed.R',
        full_dds='output/deseq2_new/data/{location}/{location}_full_dds.rds',
        samples='scripts/htseqcount_samples_full.csv'
    output:
        expressed_genes='output/deseq2_new/expressed_genes/{location}/{location}_expressed_genes.RData'
    conda:
        '../envs/conda/bioconductor-deseq2=1.30.0.yaml'
    script:
        '../scripts/deseq2/count_percentage_genes_expressed.R'

'''
snakemake --cores 1 --use-conda \
output/expressed_genes/bristol/bristol_expressed_genes.RData
'''

rule plot_percentage_expression_all_new:
    input:
        script='scripts/plot_percentage_expression_all.R',
        full_dds_altadena='output/deseq2_new/expressed_genes/altadena/altadena_expressed_genes.RData',
        full_dds_bristol='output/deseq2_new/expressed_genes/bristol/bristol_expressed_genes.RData'
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

rule de_analysis_new:
    input:
        script='scripts/deseq2/de_analysis.R',
        full_dds='output/deseq2_new/data/{location}/{location}_full_dds.rds'
    output:
        dds_res='output/deseq2_new/results/{location}/{location}_{diet1}_{diet2}.rds'
    conda:
        '../envs/conda/bioconductor-deseq2=1.30.0.yaml'
    script:
        '../scripts/deseq2/de_analysis.R'

'''
snakemake --cores 1 --use-conda output/deseq2_data/bristol/bristol_as_op50.RData
'''

rule deseq2_all_count_matrix_new:
    input:
        script='scripts/deseq2_new/deseq2_all_count_matrix_new.R',
        htseq_count=expand('from_cluster/htseq_count_export_new/{location}_{diet}_{replicate}_counts.txt', \
        location = LOCATIONS, diet=DIETS, replicate=REPLICATES)
    output:
        count_matrix='output/deseq2_new/data/all/count_matrix.rds'
    conda:
        '../envs/conda/r-tidyverse=1.2.1.yaml'
    script:
        '../scripts/deseq2_new/deseq2_all_count_matrix_new.R'

'''
snakemake --cores 1 --use-conda -R \
output/deseq2_new/data/all/count_matrix.rds
'''

rule deseq2_all_new:
    input:
        script='scripts/deseq2_new/deseq2_all_new.R',
        samples='scripts/deseq2/htseqcount_samples_full.csv',
        count_matrix='output/deseq2_new/data/all/count_matrix.rds'
    output:
        vst_dds='output/deseq2_new/data/all/all_vst_dds.rds',
        full_dds='output/deseq2_new/data/all/all_full_dds.rds'
    conda:
        '../envs/conda/bioconductor-deseq2=1.30.0.yaml'
    script:
        '../scripts/deseq2_new/deseq2_all_new.R'

'''
snakemake --cores 1 --use-conda -R \
output/deseq2_new/data/all/all_vst_dds.rds
'''

rule deseq2_pca_all_new:
    ''' PCA Plot to summarise differential expression analysis '''
    input:
        script='scripts/deseq2/all_deseq2_pca.R',
        vst_dds='output/deseq2_new/data/all/all_vst_dds.rds'
    output:
        pca_plot='output/deseq2_new/plots/all/all_deseq2_pca.tiff'
    conda:
        '../envs/conda/bioconductor-deseq2=1.30.0_r-ggplot2=3.3.1.yaml'
    script:
        '../scripts/deseq2/all_deseq2_pca.R'

'''
snakemake --cores 1 --use-conda -R \
output/deseq2_new/plots/all/all_deseq2_pca.tiff
'''
