rule trinity:
    input:
        trimmed1="polyA/QC/trimmed_fastq/{location}/{diet}/{replicate}/reads1/{location}_{diet}_{replicate}_1_trimmed.fastq",
        trimmed2="polyA/QC/trimmed_fastq/{location}/{diet}/{replicate}/reads2/{location}_{diet}_{replicate}_2_trimmed.fastq"
    output:
        "polyA/reference_free/trinity/{location}/{diet}/{replicate}/trinity/Trinity.fasta"
    conda:
        "../envs/conda/Trinity=2.11.0.yaml"
    threads: 24
    shell:
        "Trinity --seqType fq --max_memory 60G --CPU {threads} \
        --left {input.trimmed1} --right {input.trimmed2} \
        --output polyA/reference_free/trinity/{wildcards.location}/{wildcards.diet}/{wildcards.replicate}/trinity/"

snakemake --use-conda --cores 24 --snakefile rules/trinity.smk \
polyA/reference_free/trinity/bristol/hb101/rep1/trinity/Trinity.fasta

# hb101 rep1 submitted via snakemake
snakemake --cluster-config snakemake_profile/slurm.json --use-conda \
--profile snakemake_profile --cores 24 --snakefile rules/trinity.smk \
polyA/reference_free/trinity/bristol/hb101/rep1/trinity/Trinity.fasta
