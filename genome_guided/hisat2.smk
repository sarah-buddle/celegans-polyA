# rule index:
#     input:
#         "genomes/bristol_genome.fa"
#     output:
#         "genomes/bristol_genome.fa"
#     conda:
#         "../envs/conda/hisat2=2.2.1.yaml"
#     shell:
#         "hisat2-build {input} {output}"

rule hisat2:
    input:
        indexes="annotations/indexes/index",
        reads1="polyA/QC/trimmed_fastq2/{sample}_1_trimmed2.fastq",
        reads2="polyA/QC/trimmed_fastq2/{sample}_2_trimmed2.fastq"
    output:
        sam="polyA/alignment/hisat2/{sample}.sam",
        junctions="polyA/alignment/hisat2/{sample}.junctions",
        log="polyA/alignment/hisat2/{sample}.log"
    conda:
        "../envs/conda/hisat2=2.2.1.yaml"
    shell:
        "hisat2 -p 8 --dta -x {input.indexes} -1 {input.reads1} -2 {input.reads2} -S {output.sam}"

rule hisat2:
     input:
         reads1="polyA/QC/trimmed_fastq2/{sample}_1_trimmed2.fastq",
         reads2="polyA/QC/trimmed_fastq2/{sample}_2_trimmed2.fastq"
     output:
         "polyA/alignment/hisat2/{sample}.sam"
     conda:
         "../envs/conda/hisat2=2.2.1.yaml"
     shell:
         "hisat2 -p 8 --dta -x annotations/indexes/index -1 {input.reads1} -2 {input.reads2} -S {output}"

hisat2 -p 8 --dta -x annotations/indexes/index \
-1 polyA/QC/trimmed_fastq2/as_rep1_1_trimmed2.fastq \
-2 polyA/QC/trimmed_fastq2/as_rep1_2_trimmed2.fastq \
-S polyA/alignment/hisat2/as_rep1.sam

snakemake --use-conda --cores 8 --snakefile rules/alignment.smk \
polyA/alignment/hisat2/as_rep3.sam

snakemake --use-conda --cores 8 --snakefile rules/alignment.smk \
polyA/alignment/hisat2/{bp_rep1,bp_rep2,bp_rep3,\
hb101_rep1,hb101_rep2,hb101_rep3,\
m9_rep1,m9_rep2,m9_rep3,\
op50_rep1,op50_rep2,op50_rep3,\
pf_rep1,pf_rep2,pf_rep3}.sam

rule sort:
    input:
        "polyA/alignment/hisat2/{sample}.sam"
    output:
        "polyA/alignment/sorted_bam/{sample}.bam"
    conda:
        "../envs/conda/samtools=1.11.yaml"
    shell:
        "samtools sort -@ 8 -o {output} {input}"

snakemake --use-conda --cores 8 --snakefile rules/sort.smk \
polyA/alignment/sorted_bam/{as_rep1,as_rep2,as_rep3,\
bp_rep1,bp_rep2,bp_rep3,\
hb101_rep1,hb101_rep2,hb101_rep3,\
m9_rep1,m9_rep2,m9_rep3,\
op50_rep1,op50_rep2,op50_rep3,\
pf_rep1,pf_rep2,pf_rep3}.bam
