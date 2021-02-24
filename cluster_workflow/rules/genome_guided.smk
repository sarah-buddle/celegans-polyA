LOCATIONS = ['altadena', 'bristol']
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

'''
snakemake --profile ../snakemake_profile \
output/polyA/genome_guided/sorted_bam/bristol/as/rep1/bristol_as_rep1.bam.bai \
output/polyA/genome_guided/sorted_bam/bristol/as/rep2/bristol_as_rep2.bam.bai \
output/polyA/genome_guided/sorted_bam/bristol/as/rep3/bristol_as_rep3.bam.bai \
output/polyA/genome_guided/sorted_bam/bristol/bp/rep1/bristol_bp_rep1.bam.bai \
output/polyA/genome_guided/sorted_bam/bristol/bp/rep2/bristol_bp_rep2.bam.bai \
output/polyA/genome_guided/sorted_bam/bristol/bp/rep3/bristol_bp_rep3.bam.bai \
output/polyA/genome_guided/sorted_bam/bristol/hb101/rep1/bristol_hb101_rep1.bam.bai \
output/polyA/genome_guided/sorted_bam/bristol/hb101/rep2/bristol_hb101_rep2.bam.bai \
output/polyA/genome_guided/sorted_bam/bristol/hb101/rep3/bristol_hb101_rep3.bam.bai \
output/polyA/genome_guided/sorted_bam/bristol/m9/rep1/bristol_m9_rep1.bam.bai \
output/polyA/genome_guided/sorted_bam/bristol/op50/rep1/bristol_op50_rep1.bam.bai \
output/polyA/genome_guided/sorted_bam/bristol/op50/rep2/bristol_op50_rep2.bam.bai \
output/polyA/genome_guided/sorted_bam/bristol/op50/rep3/bristol_op50_rep3.bam.bai \
output/polyA/genome_guided/sorted_bam/bristol/pf/rep1/bristol_pf_rep1.bam.bai \
output/polyA/genome_guided/sorted_bam/bristol/pf/rep2/bristol_pf_rep2.bam.bai \
output/polyA/genome_guided/sorted_bam/bristol/pf/rep3/bristol_pf_rep3.bam.bai \
output/polyA/genome_guided/sorted_bam/altadena/as/rep1/altadena_as_rep1.bam.bai \
output/polyA/genome_guided/sorted_bam/altadena/as/rep2/altadena_as_rep2.bam.bai \
output/polyA/genome_guided/sorted_bam/altadena/as/rep3/altadena_as_rep3.bam.bai \
output/polyA/genome_guided/sorted_bam/altadena/bp/rep1/altadena_bp_rep1.bam.bai \
output/polyA/genome_guided/sorted_bam/altadena/bp/rep2/altadena_bp_rep2.bam.bai \
output/polyA/genome_guided/sorted_bam/altadena/bp/rep3/altadena_bp_rep3.bam.bai \
output/polyA/genome_guided/sorted_bam/altadena/hb101/rep1/altadena_hb101_rep1.bam.bai \
output/polyA/genome_guided/sorted_bam/altadena/hb101/rep2/altadena_hb101_rep2.bam.bai \
output/polyA/genome_guided/sorted_bam/altadena/hb101/rep3/altadena_hb101_rep3.bam.bai \
output/polyA/genome_guided/sorted_bam/altadena/m9/rep1/altadena_m9_rep1.bam.bai \
output/polyA/genome_guided/sorted_bam/altadena/op50/rep1/altadena_op50_rep1.bam.bai \
output/polyA/genome_guided/sorted_bam/altadena/op50/rep2/altadena_op50_rep2.bam.bai \
output/polyA/genome_guided/sorted_bam/altadena/op50/rep3/altadena_op50_rep3.bam.bai \
output/polyA/genome_guided/sorted_bam/altadena/pf/rep1/altadena_pf_rep1.bam.bai \
output/polyA/genome_guided/sorted_bam/altadena/pf/rep2/altadena_pf_rep2.bam.bai \
output/polyA/genome_guided/sorted_bam/altadena/pf/rep3/altadena_pf_rep3.bam.bai
'''

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


# was orginally run with -e option
rule stringtie:
    ''' Assemble alignments into transcripts '''
    input:
        bam='output/polyA/genome_guided/sorted_bam/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.bam',
        annotation='input/annotations/liftover/{location}/{location}_annotation.gtf'
    output:
        "output/polyA/genome_guided/stringtie/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.gtf"
    threads: 6
    conda:
        '../envs/conda/stringtie=2.1.4.yaml'
    shell:
        "stringtie {input.bam} -o {output} -p {threads} -G {input.annotation}"


'''
snakemake --profile ../snakemake_profile --allowed-rules stringtie -R \
output/polyA/genome_guided/stringtie/bristol/as/rep1/bristol_as_rep1.gtf \
output/polyA/genome_guided/stringtie/bristol/bp/rep1/bristol_bp_rep1.gtf \
output/polyA/genome_guided/stringtie/bristol/bp/rep2/bristol_bp_rep2.gtf \
output/polyA/genome_guided/stringtie/bristol/bp/rep3/bristol_bp_rep3.gtf \
output/polyA/genome_guided/stringtie/bristol/hb101/rep1/bristol_hb101_rep1.gtf \
output/polyA/genome_guided/stringtie/bristol/hb101/rep2/bristol_hb101_rep2.gtf \
output/polyA/genome_guided/stringtie/bristol/hb101/rep3/bristol_hb101_rep3.gtf \
output/polyA/genome_guided/stringtie/bristol/m9/rep1/bristol_m9_rep1.gtf \
output/polyA/genome_guided/stringtie/bristol/op50/rep1/bristol_op50_rep1.gtf \
output/polyA/genome_guided/stringtie/bristol/op50/rep2/bristol_op50_rep2.gtf \
output/polyA/genome_guided/stringtie/bristol/op50/rep3/bristol_op50_rep3.gtf \
output/polyA/genome_guided/stringtie/bristol/pf/rep1/bristol_pf_rep1.gtf \
output/polyA/genome_guided/stringtie/bristol/pf/rep2/bristol_pf_rep2.gtf \
output/polyA/genome_guided/stringtie/bristol/pf/rep3/bristol_pf_rep3.gtf
'''

rule stringtie_merge:
    input:
        rep1='output/polyA/genome_guided/stringtie/{location}/{diet}/rep1/{location}_{diet}_rep1.gtf',
        rep2='output/polyA/genome_guided/stringtie/{location}/{diet}/rep2/{location}_{diet}_rep2.gtf',
        rep3='output/polyA/genome_guided/stringtie/{location}/{diet}/rep3/{location}_{diet}_rep3.gtf'
    output:
        'output/polyA/genome_guided/stringtie/{location}/{diet}/rep123/{location}_{diet}_rep123.gtf'
    conda:
        '../envs/conda/stringtie=2.1.4.yaml'
    shell:
        'stringtie --merge {input.rep1} {input.rep2} {input.rep3} -o {output}'

'''
snakemake --profile ../snakemake_profile \
output/polyA/genome_guided/stringtie/altadena/as/rep123/altadena_as_rep123.gtf \
output/polyA/genome_guided/stringtie/altadena/bp/rep123/altadena_bp_rep123.gtf \
output/polyA/genome_guided/stringtie/altadena/hb101/rep123/altadena_hb101_rep123.gtf \
output/polyA/genome_guided/stringtie/altadena/m9/rep123/altadena_m9_rep123.gtf \
output/polyA/genome_guided/stringtie/altadena/op50/rep123/altadena_op50_rep123.gtf \
output/polyA/genome_guided/stringtie/altadena/pf/rep123/altadena_pf_rep123.gtf
'''

rule merge_stringtie_all:
    input:
        expand('output/polyA/genome_guided/stringtie/{{location}}/{diet}/rep123/{{location}}_{diet}_rep123.gtf',
        diet = DIETS)
    output:
        'output/polyA/genome_guided/stringtie/{location}/all/rep123/stringtie_{location}_all_rep123.gtf'
    conda:
        '../envs/conda/stringtie=2.1.4.yaml'
    shell:
        'stringtie --merge {input} -o {output}'

'''
snakemake --profile ../snakemake_profile \
output/polyA/genome_guided/stringtie/altadena/all/rep123/altadena_all_rep123.gtf
'''

rule stringtie_export:
    input:
        'output/polyA/genome_guided/stringtie/{location}/{diet}/rep123/{location}_{diet}_rep123.gtf'
    output:
        'output/polyA/genome_guided/stringtie_export/{location}/{diet}/rep123/stringtie_{location}_{diet}_rep123.gtf'
    shell:
        'cp {input} {output}'

'''
snakemake --profile ../snakemake_profile \
output/polyA/genome_guided/stringtie_export/altadena/as/rep123/stringtie_altadena_as_rep123.gtf \
output/polyA/genome_guided/stringtie_export/altadena/bp/rep123/stringtie_altadena_bp_rep123.gtf \
output/polyA/genome_guided/stringtie_export/altadena/hb101/rep123/stringtie_altadena_hb101_rep123.gtf \
output/polyA/genome_guided/stringtie_export/altadena/m9/rep123/stringtie_altadena_m9_rep123.gtf \
output/polyA/genome_guided/stringtie_export/altadena/op50/rep123/stringtie_altadena_op50_rep123.gtf \
output/polyA/genome_guided/stringtie_export/altadena/pf/rep123/stringtie_altadena_pf_rep123.gtf
'''

'''
# export to local machine for analysis in R
scp sb2226@172.25.11.131:/mnt/home1/miska/sb2226/workflow/output/polyA/genome_guided/stringtie/bristol/as/rep2/bristol_as_rep2.gtf \
sb2226@172.25.11.131:/mnt/home1/miska/sb2226/workflow/output/polyA/genome_guided/stringtie/bristol/as/rep3/bristol_as_rep3.gtf \
sb2226@172.25.11.131:/mnt/home1/miska/sb2226/workflow/output/polyA/genome_guided/stringtie/bristol/m9/rep2/bristol_m9_rep2.gtf \
sb2226@172.25.11.131:/mnt/home1/miska/sb2226/workflow/output/polyA/genome_guided/stringtie/bristol/m9/rep3/bristol_m9_rep3.gtf .

scp sb2226@172.25.11.131:/mnt/home1/miska/sb2226/workflow/output/polyA/genome_guided/stringtie/bristol/all/rep123/stringtie_bristol_all_rep123.gtf \
sb2226@172.25.11.131:/mnt/home1/miska/sb2226/workflow/output/polyA/genome_guided/stringtie/altadena/all/rep123/stringtie_altadena_all_rep123.gtf .

'''
