''' Differential expression analysis on data obtained from liftover annotations '''

'''
deseq2_combined.smk is used for the analysis in the report
Run genome_guided.smk on cluster and export to workflow/from_cluster
'''

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
        rlog_dds='output/deseq2/data/{location}/{location}_rlog_dds.rds',
        full_dds='output/deseq2/data/{location}/{location}_full_dds.rds'
    conda:
        '../envs/conda/bioconductor-deseq2=1.30.0.yaml'
    script:
        '../scripts/deseq2/deseq2.R'

'''
snakemake --cores 1 --use-conda -R \
output/deseq2/data/bristol/bristol_rlog_dds.rds \
output/deseq2/data/altadena/altadena_rlog_dds.rds
'''

rule deseq2_PCA:
    ''' PCA Plot to summarise differential expression analysis '''
    input:
        script='scripts/deseq2/deseq2_pca.R',
        rlog_dds='output/deseq2/data/{location}/{location}_rlog_dds.rds'
    output:
        pca_plot='output/deseq2/plots/{location}/{location}_deseq2_pca.tiff'
    conda:
        '../envs/conda/bioconductor-deseq2=1.30.0_r-ggplot2=3.3.1.yaml'
    script:
        '../scripts/deseq2/deseq2_pca.R'

'''
snakemake --cores 1 --use-conda -R \
output/deseq2/plots/bristol/bristol_deseq2_pca.tiff \
output/deseq2/plots/altadena/altadena_deseq2_pca.tiff
'''

rule deseq2_no_m9:
    ''' Creates DESeq2 object without m9/starvation diet treatment '''
    input:
        script='scripts/deseq2/deseq2_no_m9.R',
        samples='scripts/deseq2/htseqcount_samples_full.csv',
        htseq_count=expand('from_cluster/htseq_count_export/{{location}}_{diet}_{replicate}_counts.txt', \
        diet=DIETS, replicate=REPLICATES)
    output:
        rlog_dds='output/deseq2/data/{location}/{location}_no_m9_rlog_dds.rds'
    conda:
        '../envs/conda/bioconductor-deseq2=1.30.0.yaml'
    script:
        '../scripts/deseq2/deseq2_no_m9.R'


rule deseq2_PCA_no_m9:
    '''PCA plot of above '''
    input:
        script='scripts/deseq2/deseq2_pca.R',
        rlog_dds='output/deseq2/data/{location}/{location}_no_m9_rlog_dds.rds'
    output:
        pca_plot='output/deseq2/plots/{location}/{location}_deseq2_pca_no_m9.tiff'
    conda:
        '../envs/conda/bioconductor-deseq2=1.30.0_r-ggplot2=3.3.1.yaml'
    script:
        '../scripts/deseq2/deseq2_pca.R'


'''
snakemake --cores 1 --use-conda \
output/deseq2/plots/bristol/bristol_deseq2_pca_no_m9.tiff \
output/deseq2/plots/altadena/altadena_deseq2_pca_no_m9.tiff
'''

rule deseq2_all_count_matrix:
    ''' Produce count matrix for both locations '''
    input:
        script='scripts/deseq2/deseq2_all_count_matrix.R',
        htseq_count=expand('from_cluster/htseq_count_export/{location}_{diet}_{replicate}_counts.txt', \
        location = LOCATIONS, diet=DIETS, replicate=REPLICATES)
    output:
        count_matrix='output/deseq2/data/all/count_matrix.rds'
    conda:
        '../envs/conda/r-tidyverse=1.2.1.yaml'
    script:
        '../scripts/deseq2/deseq2_all_count_matrix.R'

'''
snakemake --cores 1 --use-conda -R \
output/deseq2/data/all/count_matrix.rds
'''

rule deseq2_all:
    ''' Run deseq2 on data from both locations '''
    input:
        script='scripts/deseq2/deseq2_all.R',
        samples='scripts/deseq2/htseqcount_samples_full.csv',
        count_matrix='output/deseq2/data/all/count_matrix.rds'
    output:
        vst_dds='output/deseq2/data/all/all_vst_dds.rds',
        full_dds='output/deseq2/data/all/all_full_dds.rds'
    conda:
        '../envs/conda/bioconductor-deseq2=1.30.0.yaml'
    script:
        '../scripts/deseq2/deseq2_all.R'

'''
snakemake --cores 1 --use-conda -R \
output/deseq2/data/all/all_vst_dds.rds
'''

rule deseq2_pca_all:
    ''' PCA Plot to summarise differential expression analysis for both locations '''
    input:
        script='scripts/deseq2/all_deseq2_pca.R',
        vst_dds='output/deseq2/data/all/all_vst_dds.rds'
    output:
        pca_plot='output/deseq2/plots/all/all_deseq2_pca.tiff'
    conda:
        '../envs/conda/bioconductor-deseq2=1.30.0_r-ggplot2=3.3.1.yaml'
    script:
        '../scripts/deseq2/all_deseq2_pca.R'

'''
snakemake --cores 1 --use-conda -R \
output/deseq2/plots/all/all_deseq2_pca.tiff
'''

rule deseq2_percentage_genes_expressed:
    ''' Plot percentage genes expressed for single location '''
    input:
        script='scripts/deseq2/plot_percentage_genes_expressed.R',
        full_dds='output/deseq2/data/{location}/{location}_full_dds.rds',
        samples='scripts/deseq2/htseqcount_samples_full.csv'
    output:
        percentage_plot='output/deseq2/plots/{location}/{location}_percentage_plot.tiff'
    conda:
        '../envs/conda/bioconductor-deseq2=1.30.0_r-ggplot2=3.3.1.yaml'
    script:
        '../scripts/deseq2/plot_percentage_genes_expressed.R'

'''
snakemake --cores 1 --use-conda \
output/deseq2/plots/bristol/bristol_percentage_plot.tiff \
output/deseq2/plots/altadena/altadena_percentage_plot.tiff
'''

rule count_percentage_genes_expressed:
    ''' Count percentage of total genes expressed for single location '''
    input:
        script='scripts/deseq2/count_percentage_genes_expressed.R',
        full_dds='output/deseq2/data/{location}/{location}_full_dds.rds',
        samples='scripts/deseq2/htseqcount_samples_full.csv'
    output:
        expressed_genes='output/deseq2/expressed_genes/{location}/{location}_expressed_genes.rds'
    conda:
        '../envs/conda/bioconductor-deseq2=1.30.0.yaml'
    script:
        '../scripts/deseq2/count_percentage_genes_expressed.R'

'''
snakemake --cores 1 --use-conda -R \
output/deseq2/expressed_genes/bristol/bristol_expressed_genes.rds \
output/deseq2/expressed_genes/altadena/altadena_expressed_genes.rds
'''

rule plot_percentage_expression_all:
    ''' Plot percentage of total genes expressed for both locations '''
    input:
        script='scripts/deseq2/plot_percentage_expression_all.R',
        full_dds_altadena='output/deseq2/expressed_genes/altadena/altadena_expressed_genes.rds',
        full_dds_bristol='output/deseq2/expressed_genes/bristol/bristol_expressed_genes.rds'
    output:
        percentage_plot='output/deseq2/plots/all/percentage_plot_all.tiff'
    conda:
        '../envs/conda/bioconductor-deseq2=1.30.0_r-ggplot2=3.3.1.yaml'
    script:
        '../scripts/deseq2/plot_percentage_expression_all.R'

'''
snakemake --cores 1 --use-conda -R \
output/deseq2/plots/all/percentage_plot_all.tiff
'''

rule de_analysis:
    ''' Differential expression analysis for single location '''
    input:
        script='scripts/deseq2/de_analysis.R',
        full_dds='output/deseq2/data/{location}/{location}_full_dds.rds'
    output:
        dds_res='output/deseq2/results/{location}/{location}_{diet1}_{diet2}.rds'
    conda:
        '../envs/conda/bioconductor-deseq2=1.30.0.yaml'
    script:
        '../scripts/deseq2/de_analysis.R'

'''
snakemake --cores 1 --use-conda output/deseq2/results/altadena/altadena_as_pf.rds
'''

rule de_analysis_location:
    ''' Differential expression analysis for both locations '''
    input:
        script='scripts/deseq2/de_analysis_location.R',
        full_dds='output/deseq2/data/all/all_full_dds.rds'
    output:
        dds_res='output/deseq2/results_location/all/all_{location1}_{location2}.rds'
    conda:
        '../envs/conda/bioconductor-deseq2=1.30.0.yaml'
    script:
        '../scripts/deseq2/de_analysis_location.R'

'''
snakemake --cores 1 --use-conda \
output/deseq2/results_location/all/all_bristol_altadena.rds
'''

rule proportion_new_genes_de_1:
    ''' Calculate proportion of differentially expressed genes that are accessory '''
    input:
        script='scripts/deseq2/proportion_new_genes_de.R',
        dds_res=['output/deseq2/results_location/all/all_bristol_altadena.rds', \
                'output/deseq2/results/bristol/bristol_m9_op50.rds', \
                'output/deseq2/results/bristol/bristol_m9_pf.rds', \
                'output/deseq2/results/bristol/bristol_pf_op50.rds', \
                'output/deseq2/results/bristol/bristol_pf_as.rds', \
                'output/deseq2/results/bristol/bristol_as_op50.rds', \
                'output/deseq2/results/altadena/altadena_m9_op50.rds', \
                'output/deseq2/results/altadena/altadena_m9_pf.rds', \
                'output/deseq2/results/altadena/altadena_pf_op50.rds', \
                'output/deseq2/results/altadena/altadena_pf_as.rds', \
                'output/deseq2/results/altadena/altadena_as_op50.rds'],
        count_matrix='output/deseq2/data/all/count_matrix.rds'
    output:
        proportion_new_genes_plot='output/deseq2/plots/proportion_new_genes/1/proportion_new_genes.tiff'
    conda:
        '../envs/conda/r-tidyverse=1.2.1.yaml'
    script:
        '../scripts/deseq2/proportion_new_genes_de.R'

'''
snakemake --cores 1 --use-conda -R \
output/deseq2/plots/proportion_new_genes/1/proportion_new_genes.tiff
'''
