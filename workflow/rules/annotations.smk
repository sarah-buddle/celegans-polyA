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
output/annotations/granges/liftover/altadena/altadena_liftover.RData
'''

rule granges_maker:
    input:
        script='scripts/annotations/granges_maker.R',
        maker='from_cluster/maker_annotations/{annotation_type}/{location}/{diet}/{replicate}/{annotation_type}_{location}_{diet}_{replicate}.gff'
    output:
        granges='output/annotations/granges/{annotation_type}/{location}/{diet}/{replicate}/{annotation_type}_{location}_{diet}_{replicate}.RData'
    conda:
        '../envs/conda/bioconductor-rtracklayer=1.50.0.yaml'
    script:
        '../scripts/annotations/granges_maker.R'

'''
snakemake --cores 1 --use-conda \
output/annotations/granges/soap/bristol/as/rep2/soap_bristol_as_rep2.RData
'''
