rule repeatmasker:
    ''' Mask repeats in genome - needed before running Maker '''
    input:
        genome='input/genomes/{location}/{location}_genome.fasta',
        repbase='input/celrep.ref'
    output:
        'output/polyA/reference_free/repeatmasker/{location}/{location}_genome.fatsa.masked'
    conda:
        '../envs/conda/repeatmasker=4.1.1.yaml'
    shell:
        'RepeatMasker -lib {input.repbase} {input.genome} -dir output/polyA/reference_free/repeatmasker/{wildcards.location}/'

'''
snakemake --profile ../snakemake_profile \
output/polyA/reference_free/repeatmasker/altadena/altadena_genome.fasta.masked
'''
