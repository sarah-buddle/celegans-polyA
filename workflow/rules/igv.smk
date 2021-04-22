''' Prepare annotations and read alignments to be viewed in IGV '''

LOCATIONS = ['bristol']
DIETS = ['as','bp','hb101','m9','op50','pf']
REPLICATES = ['rep1','rep2','rep3']

rule remove_contig:
    ''' Remove contig from Maker annotations '''
    input:
        'from_cluster/maker_annotations/{annotation_type}/{location}/{diet}/{replicate}/{annotation_type}_{location}_{diet}_{replicate}.{extension}'
    output:
        'output/igv/gff/no_contig/{annotation_type}/{location}/{diet}/{replicate}/{annotation_type}_{location}_{diet}_{replicate}.{extension}'
    shell:
        "sed '/contig/d' {input} > {output}"

'''
snakemake --cores 1 --use-conda \
output/igv/gff/no_contig/soap/bristol/all/rep123/soap_bristol_all_rep123.gtf
'''

rule gff_sort:
    ''' Sort annotation '''
    input:
        'output/igv/gff/no_contig/{annotation_type}/{location}/{diet}/{replicate}/{annotation_type}_{location}_{diet}_{replicate}.{extension}'
    output:
        'output/igv/gff/sorted/{annotation_type}/{location}/{diet}/{replicate}/{annotation_type}_{location}_{diet}_{replicate}.{extension}'
    conda:
        '../envs/conda/gff3sort=0.1.a1a2bc9.yaml'
    shell:
        'gff3sort.pl {input} > {output}'

'''
snakemake --cores 1 --use-conda -R \
output/igv/gff/sorted/combinedall/bristol/all/rep123/combinedall_bristol_all_rep123.gtf
'''

rule gff2bed:
    ''' Convert to bed file '''
    input:
        'output/igv/gff/sorted/{annotation_type}/{location}/{diet}/{replicate}/{annotation_type}_{location}_{diet}_{replicate}.{extension}'
    output:
        'output/igv/gff/bed/{annotation_type}/{location}/{diet}/{replicate}/{annotation_type}_{location}_{diet}_{replicate}_{extension}.bed'
    conda:
        '../envs/conda/bedops=2.4.39.yaml'
    shell:
        'gff2bed < {input} > {output}'

rule index_bed:
    ''' Index bed file '''
    input:
        'output/igv/gff/bed/{annotation_type}/{location}/{diet}/{replicate}/{annotation_type}_{location}_{diet}_{replicate}_{extension}.bed'
    output:
        'output/igv/gff/bed/{annotation_type}/{location}/{diet}/{replicate}/{annotation_type}_{location}_{diet}_{replicate}_{extension}.bed.idx'
    conda:
        '../envs/conda/igvtools=2.5.3.yaml'
    shell:
        'igvtools index {input}'

'''
snakemake --cores 1 --use-conda \
output/igv/gff/bed/soap/bristol/hb101/rep123/soap_bristol_hb101_rep123.bed.idx \
output/igv/gff/bed/soap/bristol/m9/rep123/soap_bristol_m9_rep123.bed.idx \
output/igv/gff/bed/soap/bristol/op50/rep123/soap_bristol_op50_rep123.bed.idx \
output/igv/gff/bed/soap/bristol/pf/rep123/soap_bristol_pf_rep123.bed.idx \
output/igv/gff/bed/soap/altadena/as/rep123/soap_altadena_as_rep123.bed.idx \
output/igv/gff/bed/soap/altadena/bp/rep123/soap_altadena_bp_rep123.bed.idx \
output/igv/gff/bed/soap/altadena/hb101/rep123/soap_altadena_hb101_rep123.bed.idx \
output/igv/gff/bed/soap/altadena/m9/rep123/soap_altadena_m9_rep123.bed.idx \
output/igv/gff/bed/soap/altadena/op50/rep123/soap_altadena_op50_rep123.bed.idx \
output/igv/gff/bed/soap/altadena/pf/rep123/soap_altadena_pf_rep123.bed.idx
'''

rule gtf_sort_liftover:
    ''' Sort liftover annotation '''
    input:
        'from_cluster/liftover_annotations/{location}/liftover_{location}.gtf'
    output:
        'output/igv/liftover/sorted/{location}/liftover_{location}.gtf'
    conda:
        '../envs/conda/gff3sort=0.1.a1a2bc9.yaml'
    shell:
        'gff3sort.pl {input} > {output}'

'''
snakemake --cores 1 --use-conda \
output/igv/liftover/sorted/altadena/liftover_altadena.gtf
'''

rule gff2bed_liftover:
    ''' Convert to bed file '''
    input:
        'output/igv/liftover/sorted/{location}/liftover_{location}.gtf'
    output:
        'output/igv/liftover/bed/{location}/liftover_{location}.bed'
    conda:
        '../envs/conda/bedops=2.4.39.yaml'
    shell:
        'gff2bed < {input} > {output}'

rule index_bed_liftover:
    ''' Index bed file '''
    input:
        'output/igv/liftover/bed/{location}/liftover_{location}.bed'
    output:
        'output/igv/liftover/bed/{location}/liftover_{location}.bed.idx'
    conda:
        '../envs/conda/igvtools=2.5.3.yaml'
    shell:
        'igvtools index {input}'

'''
snakemake --cores 1 --use-conda \
output/igv/liftover/bed/bristol/liftover_bristol.bed.idx
'''

rule bam_lists:
    ''' Create lists of bam files to combine biological replicates into single track '''
    output:
        'output/igv/bam_lists/{location}/{diet}/{location}_{diet}.bam.list'
    shell:
        'echo /Users/Sarah/Documents/sorted_bam/{wildcards.location}/{wildcards.diet}/rep1/{wildcards.location}_{wildcards.diet}_rep1.bam >> {output};'
        'echo /Users/Sarah/Documents/sorted_bam/{wildcards.location}/{wildcards.diet}/rep2/{wildcards.location}_{wildcards.diet}_rep2.bam >> {output};'
        'echo /Users/Sarah/Documents/sorted_bam/{wildcards.location}/{wildcards.diet}/rep3/{wildcards.location}_{wildcards.diet}_rep3.bam >> {output}'

'''
snakemake --cores 1 --use-conda \
output/igv/bam_lists/bristol/as/bristol_as.bam.list \
output/igv/bam_lists/bristol/bp/bristol_bp.bam.list \
output/igv/bam_lists/bristol/hb101/bristol_hb101.bam.list \
output/igv/bam_lists/bristol/m9/bristol_m9.bam.list \
output/igv/bam_lists/bristol/op50/bristol_op50.bam.list \
output/igv/bam_lists/bristol/pf/bristol_pf.bam.list \
output/igv/bam_lists/altadena/as/altadena_as.bam.list \
output/igv/bam_lists/altadena/bp/altadena_bp.bam.list \
output/igv/bam_lists/altadena/hb101/altadena_hb101.bam.list \
output/igv/bam_lists/altadena/m9/altadena_m9.bam.list \
output/igv/bam_lists/altadena/op50/altadena_op50.bam.list \
output/igv/bam_lists/altadena/pf/altadena_pf.bam.list
'''
