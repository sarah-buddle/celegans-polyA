rule extract_ss_exon:
    ''' Extract splice sites and exons from the annotation for use in indexing '''
    input:
        'input/annotations/{location}/{location}_annotation.gtf'
    output:
        annotation='polyA/genome_guided/annotations/{location}/{location}_annotation',
        ss='polyA/genome_guided/annotations/{location}/{location}.ss',
        exon='polyA/genome_guided/annotations/{location}/{location}.exon'
    conda:
        'envs/conda/hisat2=2.2.1.yaml'
    shell:
        'cp {input} polyA/genome_guided/annotations/{location}/;',
        'extract_splice_sites.py {input} > {output.ss};',
        'extract_exons.py {input} > {output.exon}'

rule index_genome:
    ''' Index genome ready for hisat2 '''
    input:
        ss='polyA/genome_guided/annotations/{location}/{location}.ss',
        exon='polyA/genome_guided/annotations/{location}/{location}.exon',
        genome='input/genome/{location}_genome.fa'
    output:
        'polyA/genome_guided/annotations/{location}/indexes'
    conda:
        'envs/conda/hisat2=2.2.1.yaml'
    shell:
        'hisat2-build --ss {input.ss} --exon {input.exon} {input.genome} {output}'

rule hisat2:
    ''' Align reads to the genome '''
    input:
        indexes='polyA/genome_guided/annotations/{location}/indexes',
        trimmed1='polyA/QC/trimmed_fastq/{location}/{diet}/{replicate}/reads1/{location}_{diet}_{replicate}_1_trimmed.fastq',
        trimmed2='polyA/QC/trimmed_fastq/{location}/{diet}/{replicate}/reads2/{location}_{diet}_{replicate}_2_trimmed.fastq'
    output:
        'polyA/genome_guided/hisat2/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.sam'
    conda:
        'envs/conda/hisat2=2.2.1.yaml'
    shell:
        'hisat2 -p 8 --dta -x {input.indexes} -1 {input.trimmed1} -2 {input.trimmed2} -S {output}'

rule sam_to_sorted_bam:
    ''' Convert sam file to bam and sort '''
    input:
        'polyA/genome_guided/hisat2/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.sam'
    output:
        'polyA/genome_guided/sorted_bam/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.bam'
    conda:
        'envs/conda/samtools=1.11.yaml'
    shell:
        "samtools sort -@ 8 -o {output} {input}"

rule stringtie:
    ''' Assemble alignments into transcripts '''
    input:
        bam='polyA/genome_guided/sorted_bam/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.bam'
        annotation='input/annotations/{location}/{location}_annotation.gtf'
    output:
        "polyA/genome_guided/stringtie/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.gtf"
    conda:
        'envs/conda/stringtie=2.1.4.yaml'
    shell:
        "stringtie {input.bam} -o {output} -p 8 -G {input.annotation} -e"

rule make_sample_list:
    ''' Make text file of sample names and stringtie file paths for prepDE.py script'''
    output:
        sample_paths='polyA/genome_guided/deseq2_prep/sample_paths.txt',
        sample_names='polyA/genome_guided/deseq2_prep/sample_names.txt',
        sample_list='polyA/genome_guided/deseq2_prep/sample_list.txt'
    shell:
        'find polyA/genome_guided/stringtie -type f | sort > {output.sample_paths};'
        'find polyA/genome_guided/stringtie -type f -exec basename {{}} \; | sort | rev | cut -c5- | rev > {output.sample_names};'
        'paste {output.sample_names} {output.sample_paths} > {output.sample_list}'

'''
snakemake --cores 1 --snakefile rules/make_sample_list.smk polyA/genome_guided/deseq2_prep/sample_list.txt \
polyA/genome_guided/deseq2_prep/sample_names.txt polyA/genome_guided/deseq2_prep/sample_paths.txt
'''

rule prepDE:
    ''' Produce gene and transcript count matrices for DESeq2 '''
    input:
        'polyA/genome_guided/deseq2_prep/sample_list.txt'
    output:
        gene_count='polyA/genome_guided/deseq2_prep/gene_count_matrix.csv',
        transcript_count='polyA/genome_guided/deseq2_prep/transcript_count_matrix.csv'
    shell:
        'python2.7 scripts/prepDE.py -i {input} -g {output.gene_count} -t {output.transcript_count}'

'''
snakemake --cores 1 --snakefile rules/prepDE.smk \
polyA/genome_guided/deseq2_prep/gene_count_matrix.csv \
polyA/genome_guided/deseq2_prep/transcript_count_matrix.csv
'''

'''
# export to local machine for analysis in R
exit
cd
scp sb2226@172.25.11.131:/mnt/home1/miska/sb2226/polyA/genome_guided/deseq2_prep/gene_count_matrix.csv \
/Users/Sarah/OneDrive/Documents/Uni/III/Project/Scripts/polyA/genome_guided

scp sb2226@172.25.11.131:/mnt/home1/miska/sb2226/polyA/genome_guided/deseq2_prep/transcript_count_matrix.csv \
/Users/Sarah/OneDrive/Documents/Uni/III/Project/Scripts/polyA/genome_guided
'''
