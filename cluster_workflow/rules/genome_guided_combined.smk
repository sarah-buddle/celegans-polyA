''' Counts read mapping to each gene in new combined annotation '''

LOCATIONS = ['altadena', 'bristol']
DIETS = ['as','bp','hb101','m9','op50','pf']
REPLICATES = ['rep1','rep2','rep3']
INDEXES = ['1','2','3','4','5','6','7','8']

rule htseq_count_combined:
    ''' Counts read mapping to each gene in new combined annotation '''
    input:
        bam='output/polyA/genome_guided/sorted_bam/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.bam',
        annotation='input/annotations/combinedall/{location}/combinedall_{location}.gtf'
    output:
        'output/polyA/genome_guided_combined/htseq_count/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}_counts.txt'
    conda:
        '../envs/conda/htseq=0.12.4.yaml'
    shell:
        'htseq-count --format=bam --stranded=no --type=gene {input.bam} {input.annotation} > {output}'

'''
snakemake --profile ../snakemake_profile -R \
output/polyA/genome_guided_combined/htseq_count/altadena/as/rep1/altadena_as_rep1_counts.txt
'''

rule htseq_count_export_combined:
    ''' Moves count data to single file for export to local machine '''
    input:
        expand('output/polyA/genome_guided_combined/htseq_count/{{location}}/{diet}/{replicate}/{{location}}_{diet}_{replicate}_counts.txt', \
        diet=DIETS, replicate=REPLICATES)
    output:
        'output/polyA/genome_guided_combined/htseq_count/htseq_count_export/{location}_{diet}_{replicate}_counts.txt'
    shell:
        'cp {input} output/polyA/genome_guided_combined/htseq_count/htseq_count_export/'

'''
snakemake --profile ../snakemake_profile \
output/polyA/genome_guided_combined/htseq_count/htseq_count_export/bristol_{as_rep2,as_rep3,\
bp_rep1,bp_rep2,bp_rep3,\
hb101_rep1,hb101_rep2,hb101_rep3,\
m9_rep1,m9_rep2,m9_rep3,\
op50_rep1,op50_rep2,op50_rep3,\
pf_rep1,pf_rep2,pf_rep3}_counts.txt \
output/polyA/genome_guided_combined/htseq_count/htseq_count_export/bristol_{bp_rep1,bp_rep2,bp_rep3,\
hb101_rep1,hb101_rep2,hb101_rep3,\
m9_rep1,m9_rep2,m9_rep3,\
op50_rep1,op50_rep2,op50_rep3,\
pf_rep1,pf_rep2,pf_rep3}_counts.txt
'''

'''
Export output/polyA/genome_guided_combined/htseq_count/htseq_count_export/ to
workflow/from_cluster/htseq_count_export_combined on local machine
'''
