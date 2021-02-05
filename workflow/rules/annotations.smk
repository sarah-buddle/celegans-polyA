rule granges_liftover:
    input:
        script='scripts/annotations/granges_liftover.R',
        liftover='from_cluster/liftover_annotations/{location}/liftover_{location}.gtf'
    output:
        granges='output/annotations/granges/liftover/{location}/{location}_liftover.RData'
    conda:
        '../envs/conda/bioconductor-rtracklayer=1.50.0.yaml'
    script:
        '../scripts/annotations/granges_liftover.R'

'''
snakemake --cores 1 --use-conda \
output/annotations/granges/liftover/bristol/bristol_liftover.RData
'''
