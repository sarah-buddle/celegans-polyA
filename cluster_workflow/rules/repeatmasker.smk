# the moving files part at the end doesn't work
rule repeatmasker:
    input:
        genome='input/genomes/{location}/{location}_genome.fasta',
        repbase='input/celrep.ref'
    output:
        'polyA/reference_free/repeatmasker/{location}/{location}_genome_masked.fasta'
    conda:
        '../envs/conda/repeatmasker=4.1.1.yaml'
    shell:
        'RepeatMasker -lib {input.repbase} {input.genome};'
        'mv input/genomes/{wildcards.location}/{wildcards.location}_genome.fasta.masked polyA/reference_free/repeatmasker/{wildcards.location}/;'
        'mv polyA/reference_free/repeatmasker/{wildcards.location}/{wildcards.location}_genome.fasta.masked {output}'

snakemake --cluster-config scripts/slurm.json --profile snakemake_profile \
--use-conda --cores 24 --snakefile rules/repeatmasker.smk \
polyA/reference_free/repeatmasker/bristol/bristol_genome_masked.fasta

snakemake --cores 1 --use-conda --snakefile rules/repeatmasker.smk \
polyA/reference_free/repeatmasker/bristol/bristol_genome_masked.fasta

# this command used to move instead
mv input/genomes/bristol/bristol_genome.fasta.masked \
input/genomes/bristol/bristol_genome.fasta.cat.gz \
input/genomes/bristol/bristol_genome.fasta.ori.out \
input/genomes/bristol/bristol_genome.fasta.out \
input/genomes/bristol/bristol_genome.fasta.tbl \
polyA/reference_free/repeatmasker/bristol
