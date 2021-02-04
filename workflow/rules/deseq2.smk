LOCATIONS = ['bristol']
DIETS = ['as','bp','hb101','m9','op50','pf']
REPLICATES = ['rep1','rep2','rep3']

rule deseq2:
    input:
        expand('from_cluster/htseq_count_export/{location}_{diet}_{replicate}_counts.txt', \
        location=LOCATIONS, diet=DIETS, replicate=REPLICATES)
    output:
        'output/deseq2/pca_plot.png'
    conda:
        '../envs/conda/bioconductor-deseq2=1.30.0.yaml'
    script:
        '../scripts/deseq2.R'

snakemake --use-conda --cores 1 --snakefile rules/deseq2.smk \
output/deseq2/pca_plot.png
