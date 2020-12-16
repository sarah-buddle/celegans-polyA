rule htseq_count:
    ''' counts reads mapping to each gene in reference annotation '''
    input:
        bam='polyA/genome_guided/sorted_bam/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.bam',
        annotation='input/annotations/{location}/{location}_annotation.gtf'
    output:
        'polyA/genome_guided/htseq_count/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}_counts.txt'
    conda:
        '../envs/conda/htseq=0.12.4.yaml'
    shell:
        'htseq-count {input.bam} {input.annotation} > {output}'

snakemake --use-conda --cores 1 --snakefile rules/htseq_count.smk \
polyA/genome_guided/htseq_count/bristol/as/rep1/bristol_as_rep1_counts.txt

# not tested yet - see htseq_count.sh for script used
snakemake --jobs 18 --cluster-config --cluster "sbatch -n 24 -N 1 -p IACT --job-name htseq_count_snakemake \
--mail-user=sb2226@cam.ac.uk" \
--use-conda --snakefile rules/htseq_count.smk \
polyA/genome_guided/htseq_count/bristol/{as,bp,hb101,m9,op50,pf}/rep{1,2,3}/bristol_{as,bp,hb101,m9,op50,pf}_rep{1,2,3}_counts.txt

rule htseq_count_export:
    ''' moves count data to single file for export to local machine '''
    input:
        'polyA/genome_guided/htseq_count/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}_counts.txt'
    output:
        'polyA/genome_guided/htseq_count/htseq_count_export/{location}_{diet}_{replicate}_counts.txt'
    shell:
        'cp {input} polyA/genome_guided/htseq_count/htseq_count_export/'

snakemake --cores 1 --snakefile rules/htseq_count_export \
polyA/genome_guided/htseq_count/htseq_count_export/bristol_{as_rep1,as_rep2,as_rep3,\
bp_rep1,bp_rep2,bp_rep3,\
hb101_rep1,hb101_rep2,hb101_rep3,\
m9_rep1,m9_rep2,m9_rep3,\
op50_rep1,op50_rep2,op50_rep3,\
pf_rep1,pf_rep2,pf_rep3}_counts.txt

# export to local machine
scp -r sb2226@172.25.11.131:/mnt/home1/miska/sb2226/polyA/genome_guided/htseq_count/htseq_count_export \
/Users/Sarah/OneDrive/Documents/Uni/III/Project/Scripts/polyA/genome_guided
