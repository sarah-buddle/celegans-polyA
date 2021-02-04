INDEXES = ['1','2','3','4','5','6','7','8']

rule extract_unmapped_reads:
    ''' Extract unmapped reads from alignments '''
    input:
        'output/polyA/genome_guided/sorted_bam/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.bam'
    output:
        'output/polyA/decontamination/unmapped_bam/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.bam'
    conda:
        '../envs/conda/samtools=1.11.yaml'
    shell:
        'samtools view -b -f 4 {input} > {output}'

'''
snakemake --profile ../snakemake_profile \
output/polyA/decontamination/unmapped_bam/bristol/as/rep3/bristol_as_rep3.bam
'''

rule unmapped_bam_to_fastq:
    ''' Convert bam files to fastq '''
    input:
        'output/polyA/decontamination/unmapped_bam/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.bam'
    output:
        'output/polyA/decontamination/unmapped_unsorted_fastq/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.fastq'
    conda:
        '../envs/conda/samtools=1.11.yaml'
    shell:
        'samtools bam2fq {input} > {output}'

rule unmapped_sort_fastq:
    ''' Sort fastq files '''
    input:
        'output/polyA/decontamination/unmapped_unsorted_fastq/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.fastq'
    output:
        'output/polyA/decontamination/unmapped_sorted_fastq/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.fastq'
    conda:
        '../envs/conda/fastq-tools=0.8.3.yaml'
    shell:
        'fastq-sort {input} > {output}'

'''
snakemake --cluster-config ../snakemake_profile/slurm.json --use-conda \
--profile ../snakemake_profile --cores 8 \
output/polyA/decontamination/unmapped_sorted_fastq/bristol/as/rep1/bristol_as_rep1.fastq
'''

rule velveth:
    ''' De novo assembly of unmapped reads '''
    input:
        'output/polyA/decontamination/unmapped_sorted_fastq/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.fastq'
    output:
        'output/polyA/decontamination/velvet/{location}/{diet}/{replicate}/Roadmaps',
        'output/polyA/decontamination/velvet/{location}/{diet}/{replicate}/Sequences',
        'output/polyA/decontamination/velvet/{location}/{diet}/{replicate}/Log'
    conda:
        '../envs/conda/velvet=1.2.10.yaml'
    shell:
        '../anaconda3/envs/velvet=1.2.10/bin/velveth \
        output/polyA/decontamination/velvet/{wildcards.location}/{wildcards.diet}/{wildcards.replicate}/ \
        21 -fastq -shortPaired {input}'

'''
snakemake --cluster-config ../snakemake_profile/slurm.json --use-conda \
--profile ../snakemake_profile --cores 8 \
output/polyA/decontamination/velvet/bristol/as/rep1/bristol_as_rep1.Roadmaps
'''

rule velvetg:
    ''' De novo assembly of unmapped reads '''
    input:
        'output/polyA/decontamination/velvet/{location}/{diet}/{replicate}/Roadmaps',
        'output/polyA/decontamination/velvet/{location}/{diet}/{replicate}/Sequences',
        'output/polyA/decontamination/velvet/{location}/{diet}/{replicate}/Log'
    output:
        'output/polyA/decontamination/velvet/{location}/{diet}/{replicate}/contigs.fa'
    conda:
        '../envs/conda/velvet=1.2.10.yaml'
    shell:
        '../anaconda3/envs/velvet=1.2.10/bin/velvetg \
        output/polyA/decontamination/velvet/{wildcards.location}/{wildcards.diet}/{wildcards.replicate}/ \
        -cov_cutoff 4'

'''
snakemake --cluster-config ../snakemake_profile/slurm.json --use-conda \
--profile ../snakemake_profile --cores 8 \
output/polyA/decontamination/velvet/bristol/as/rep1/contigs.fa
'''

rule index_velvet:
    ''' Index velvet output ready for hisat2 '''
    input:
        'output/polyA/decontamination/velvet/{location}/{diet}/{replicate}/contigs.fa'
    output:
        expand('output/polyA/decontamination/velvet/{{location}}/{{diet}}/{{replicate}}/indexes/indexes.{index}.ht2', \
        index=INDEXES),
    conda:
        '../envs/conda/hisat2=2.2.1.yaml'
    shell:
        'hisat2-build {input} output/polyA/decontamination/velvet/{wildcards.location}/{wildcards.diet}/{wildcards.replicate}/indexes;'
        'cd output/polyA/decontamination/velvet/{wildcards.location}/{wildcards.diet}/{wildcards.replicate}/;'
        'mkdir indexes;'
        'mv indexes.*.ht2 indexes'

'''
snakemake --cluster-config ../snakemake_profile/slurm.json --use-conda \
--profile ../snakemake_profile --cores 2 \
output/polyA/decontamination/velvet/bristol/as/rep1/indexes
'''

rule hisat2_velvet:
    ''' Align reads to the genome '''
    input:
        indexes=expand('output/polyA/decontamination/velvet/{{location}}/{{diet}}/{{replicate}}/indexes/indexes.{index}.ht2', \
        index=INDEXES),
        trimmed1='output/polyA/QC/trimmed_fastq/{location}/{diet}/{replicate}/reads1/{location}_{diet}_{replicate}_1_trimmed.fastq',
        trimmed2='output/polyA/QC/trimmed_fastq/{location}/{diet}/{replicate}/reads2/{location}_{diet}_{replicate}_2_trimmed.fastq'
    output:
        'output/polyA/decontamination/hisat2_velvet/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.sam'
    conda:
        '../envs/conda/hisat2=2.2.1.yaml'
    threads: 8
    shell:
        'hisat2 --dta -x output/polyA/decontamination/velvet/{wildcards.location}/{wildcards.diet}/{wildcards.replicate}/indexes/indexes \
        -1 {input.trimmed1} -2 {input.trimmed2} -S {output}'

'''
snakemake --profile ../snakemake_profile \
output/polyA/decontamination/hisat2_velvet/bristol/as/rep1/bristol_as_rep1.sam
'''

rule sam_to_sorted_bam_2:
    ''' Convert sam file produced by hisat2 from velvet to bam and sort '''
    input:
        'output/polyA/decontamination/hisat2_velvet/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.sam'
    output:
        'output/polyA/decontamination/sorted_bam_2/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.bam'
    conda:
        '../envs/conda/samtools=1.11.yaml'
    shell:
        "samtools sort -@ 8 -o {output} {input}"

rule extract_unmapped_reads_2:
    ''' Extract unmapped reads from alignments '''
    input:
        'output/polyA/decontamination/sorted_bam_2/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.bam'
    output:
        'output/polyA/decontamination/unmapped_bam_2/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.bam'
    conda:
        '../envs/conda/samtools=1.11.yaml'
    shell:
        'samtools view -b -f 4 {input} > {output}'

'''
snakemake --profile ../snakemake_profile \
output/polyA/decontamination/unmapped_bam/bristol/as/rep3/bristol_as_rep3.bam
'''

rule unmapped_bam_to_fastq_2:
    ''' Convert bam files to fastq '''
    input:
        'output/polyA/decontamination/unmapped_bam_2/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.bam'
    output:
        'output/polyA/decontamination/unmapped_unsorted_fastq_2/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.fastq'
    conda:
        '../envs/conda/samtools=1.11.yaml'
    shell:
        'samtools bam2fq {input} > {output}'

rule unmapped_sort_fastq_2:
    ''' Sort fastq files '''
    input:
        'output/polyA/decontamination/unmapped_unsorted_fastq_2/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.fastq'
    output:
        'output/polyA/decontamination/unmapped_sorted_fastq_2/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.fastq'
    conda:
        '../envs/conda/fastq-tools=0.8.3.yaml'
    shell:
        'fastq-sort {input} > {output}'

rule deinterleave_decontamination:
    ''' Separate paired end reads '''
    input:
        'output/polyA/decontamination/unmapped_sorted_fastq_2/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.fastq'
    output:
        reads1 = 'output/polyA/decontamination/deinterleaved_fastq/{location}/{diet}/{replicate}/reads1/{location}_{diet}_{replicate}_1.fastq',
        reads2 = 'output/polyA/decontamination/deinterleaved_fastq/{location}/{diet}/{replicate}/reads2/{location}_{diet}_{replicate}_2.fastq'
    conda:
        '../envs/conda/seqtk=1.3.yaml'
    shell:
        'seqtk seq -1 {input} > {output.reads1}; seqtk seq -2 {input} > {output.reads2}'

'''
snakemake --profile ../snakemake_profile \
output/polyA/decontamination/deinterleaved_fastq/bristol/as/rep1/reads1/bristol_as_rep1_1.fastq \
output/polyA/decontamination/deinterleaved_fastq/bristol/as/rep1/reads2/bristol_as_rep1_2.fastq
'''
