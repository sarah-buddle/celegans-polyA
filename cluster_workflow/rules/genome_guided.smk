LOCATIONS = ['altadena']
DIETS = ['as','bp','hb101','m9','op50','pf']
REPLICATES = ['rep1','rep2','rep3']
INDEXES = ['1','2','3','4','5','6','7','8']

rule extract_ss_exon:
    ''' Extract splice sites and exons from the annotation for use in indexing '''
    input:
        'input/annotations/liftover/{location}/{location}_annotation.gtf'
    output:
        annotation='output/polyA/genome_guided/annotations/{location}/{location}_annotation.gtf',
        ss='output/polyA/genome_guided/annotations/{location}/{location}.ss',
        exon='output/polyA/genome_guided/annotations/{location}/{location}.exon'
    conda:
        '../envs/conda/hisat2=2.2.1.yaml'
    shell:
        #'cp {input} output/polyA/genome_guided/annotations/{wildcards.location};'
        'extract_splice_sites.py {input} > {output.ss};'
        'extract_exons.py {input} > {output.exon}'

'''
snakemake --cluster-config snakemake_profile/slurm.json --use-conda \
--profile snakemake_profile --cores 2 --snakefile rules/genome_guided_all.smk \
output/polyA/genome_guided/annotations/altadena/altadena.{ss,exon}

'''

rule index_genome:
    ''' Index genome ready for hisat2 '''
    input:
        ss='output/polyA/genome_guided/annotations/{location}/{location}.ss',
        exon='output/polyA/genome_guided/annotations/{location}/{location}.exon',
        genome='input/genomes/{location}/{location}_genome.fasta'
    output:
        'output/polyA/genome_guided/annotations/{location}/indexes'
    conda:
        '../envs/conda/hisat2=2.2.1.yaml'
    shell:
        'hisat2-build --ss {input.ss} --exon {input.exon} {input.genome} {output};'
        'mkdir indexes;'
        'mv indexes.*.ht2 indexes'

'''
snakemake --cluster-config snakemake_profile/slurm.json --use-conda \
--profile snakemake_profile --cores 2 --snakefile rules/genome_guided_all.smk \
output/polyA/genome_guided/annotations/altadena/indexes
'''

rule hisat2:
    ''' Align reads to the genome '''
    input:
        indexes=expand('output/polyA/genome_guided/annotations/{{location}}/indexes/indexes.{index}.ht2', \
        index=INDEXES),
        trimmed1='output/polyA/QC/trimmed_fastq/{location}/{diet}/{replicate}/reads1/{location}_{diet}_{replicate}_1_trimmed.fastq',
        trimmed2='output/polyA/QC/trimmed_fastq/{location}/{diet}/{replicate}/reads2/{location}_{diet}_{replicate}_2_trimmed.fastq'
    output:
        'output/polyA/genome_guided/hisat2/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.sam'
    conda:
        '../envs/conda/hisat2=2.2.1.yaml'
    shell:
        'hisat2 --dta -x output/polyA/genome_guided/annotations/{wildcards.location}/indexes/indexes -1 {input.trimmed1} -2 {input.trimmed2} -S {output}'

'''
snakemake --cluster-config snakemake_profile/slurm.json --use-conda \
--profile snakemake_profile --cores 8 --snakefile rules/genome_guided_all.smk \
output/polyA/genome_guided/hisat2/altadena/as/rep1/altadena_as_rep1.sam
'''

rule sam_to_sorted_bam:
    ''' Convert sam file to bam and sort '''
    input:
        'output/polyA/genome_guided/hisat2/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.sam'
    output:
        'output/polyA/genome_guided/sorted_bam/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.bam'
    conda:
        '../envs/conda/samtools=1.11.yaml'
    shell:
        "samtools sort -@ 8 -o {output} {input}"

'''
snakemake --cluster-config snakemake_profile/slurm.json --use-conda \
--profile snakemake_profile --cores 8 --snakefile rules/genome_guided_all.smk \
output/polyA/genome_guided/sorted_bam/altadena/as/rep1/altadena_as_rep1.bam
'''

rule index_bam:
    input:
        'output/polyA/genome_guided/sorted_bam/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.bam'
    output:
        'output/polyA/genome_guided/sorted_bam/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.bam.bai'
    conda:
        '../envs/conda/samtools=1.11.yaml'
    shell:
        'samtools index {input}'

rule htseq_count:
    ''' counts reads mapping to each gene in reference annotation '''
    input:
        bam='output/polyA/genome_guided/sorted_bam/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.bam',
        annotation='output/polyA/genome_guided/annotations/{location}/{location}_annotation.gtf'
    output:
        'output/polyA/genome_guided/htseq_count/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}_counts.txt'
    conda:
        '../envs/conda/htseq=0.12.4.yaml'
    shell:
        'htseq-count {input.bam} {input.annotation} > {output}'

'''
snakemake --cluster-config snakemake_profile/slurm.json --use-conda \
--profile snakemake_profile --cores 8 --snakefile rules/genome_guided_all.smk \
output/polyA/genome_guided/htseq_count/altadena/as/rep1/altadena_as_rep1_counts.txt
'''

rule htseq_count_export:
    ''' moves count data to single file for export to local machine '''
    input:
        expand('output/polyA/genome_guided/htseq_count/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}_counts.txt', \
        location=LOCATIONS, diet=DIETS, replicate=REPLICATES)
    output:
        'output/polyA/genome_guided/htseq_count/htseq_count_export/{location}_{diet}_{replicate}_counts.txt'
    shell:
        'cp {input} output/polyA/genome_guided/htseq_count/htseq_count_export/'

'''
snakemake --cluster-config snakemake_profile/slurm.json --use-conda \
--profile snakemake_profile --cores 24 --snakefile rules/genome_guided_all.smk \
output/polyA/genome_guided/htseq_count/htseq_count_export/altadena_{bp_rep1,bp_rep2,bp_rep3,\
hb101_rep1,hb101_rep2,hb101_rep3,\
m9_rep1,m9_rep2,m9_rep3,\
op50_rep1,op50_rep2,op50_rep3,\
pf_rep1,pf_rep2,pf_rep3}_counts.txt
'''

'''
scp -r sb2226@172.25.11.131:/mnt/home1/miska/sb2226/output/polyA/genome_guided/htseq_count/htseq_count_export
/Users/Sarah/OneDrive/Documents/Uni/III/Project/from_cluster/htseq_count_export
'''


rule stringtie:
    ''' Assemble alignments into transcripts '''
    input:
        bam='output/polyA/genome_guided/sorted_bam/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.bam',
        annotation='output/polyA/genome_guided/annotations/{location}/{location}_annotation.gtf'
    output:
        "output/polyA/genome_guided/stringtie/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.gtf"
    conda:
        '../envs/conda/stringtie=2.1.4.yaml'
    shell:
        "stringtie {input.bam} -o {output} -p 8 -G {input.annotation} -e"


'''
snakemake --cluster-config snakemake_profile/slurm.json --use-conda \
--profile snakemake_profile --cores 24 --snakefile rules/genome_guided_all.smk \
output/polyA/genome_guided/htseq_count/htseq_count_export/altadena_as_rep1_counts.txt
'''

'''
# export to local machine for analysis in R
exit
cd
scp sb2226@172.25.11.131:/mnt/home1/miska/sb2226/output/polyA/genome_guided/deseq2_prep/gene_count_matrix.csv \
/Users/Sarah/OneDrive/Documents/Uni/III/Project/Scripts/output/polyA/genome_guided

scp sb2226@172.25.11.131:/mnt/home1/miska/sb2226/output/polyA/genome_guided/deseq2_prep/transcript_count_matrix.csv \
/Users/Sarah/OneDrive/Documents/Uni/III/Project/Scripts/output/polyA/genome_guided
'''
