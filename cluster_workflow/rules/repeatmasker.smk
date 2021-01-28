# the moving files part at the end doesn't work
rule repeatmasker:
    input:
        genome='input/genomes/{location}/{location}_genome.fasta',
        repbase='input/celrep.ref'
    output:
        'output/polyA/reference_free/repeatmasker/{location}/{location}_genome_masked.fasta'
    conda:
        '../envs/conda/repeatmasker=4.1.1.yaml'
    threads: 8
    shell:
        'RepeatMasker -lib {input.repbase} {input.genome} -dir output/polyA/reference_free/repeatmasker/{wildcards.location}/'

'''
snakemake --cluster-config ../snakemake_profile/slurm.json \
--profile ../snakemake_profile --use-conda --cores 8 \
output/polyA/reference_free/repeatmasker/altadena/altadena_genome.fasta.masked
'''
