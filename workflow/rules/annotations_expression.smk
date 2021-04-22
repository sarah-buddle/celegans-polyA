''' Measure expression of new genes '''

''' Run deseq2_combined.smk '''

rule new_gene_expression:
    ''' Use data from deseq2 to measure expression of new genes '''
    input:
        script='scripts/annotations/new_gene_expression.R',
        full_dds='output/deseq2_combined/data/{location}/{location}_full_dds.rds',
        allnewgenes='output/annotations/combined_annotations/combinedall/{location}/allnewgenes_{location}.rds'
    output:
        new_gene_counts='output/annotations/new_gene_expression/data/{location}/new_gene_counts_{location}.rds',
        expression_counts='output/annotations/new_gene_expression/data/{location}/new_gene_expression_{location}.rds'
    conda:
        '../envs/conda/bioconductor-deseq2=1.30.0_bioconductor-plyranges=1.10.0.yaml'
    script:
        '../scripts/annotations/new_gene_expression.R'

'''
snakemake --cores 1 --use-conda -R \
output/annotations/new_gene_expression/data/altadena/new_gene_expression_altadena.rds \
output/annotations/new_gene_expression/data/bristol/new_gene_expression_bristol.rds
'''

rule plot_new_gene_expression:
    ''' Plot of annotation type and expression of new genes '''
    input:
        script='scripts/annotations/plot_new_gene_expression.R',
        expression_counts_bristol='output/annotations/new_gene_expression/data/bristol/new_gene_expression_bristol.rds',
        expression_counts_altadena='output/annotations/new_gene_expression/data/altadena/new_gene_expression_altadena.rds'
    output:
        new_gene_expression_plot='output/annotations/new_gene_expression/plots/new_gene_expression_plot.tiff'
    conda:
        '../envs/conda/r-ggplot2=3.3.1.yaml'
    script:
        '../scripts/annotations/plot_new_gene_expression.R'

'''
snakemake --cores 1 --use-conda -R --allowed-rules plot_new_gene_expression \
output/annotations/new_gene_expression/plots/new_gene_expression_plot.tiff
'''
