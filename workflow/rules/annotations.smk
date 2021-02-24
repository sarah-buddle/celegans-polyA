LOCATIONS = ['bristol']
DIETS = ['as','bp','hb101','m9','op50','pf']
REPLICATES = ['rep1','rep2','rep3']
CHROMOSOMES = ['chrI','chrII','chrIII','chrIV','chrV','chrX']

'''
scp -r sb2226@172.25.11.131://mnt/home1/miska/sb2226/workflow/output/polyA/reference_free/maker_soap_export \
from_cluster/maker_annotations/soap
'''

# rule rename_maker:
#     input:
#         'from_cluster/maker_annotations/soap/{location}/{diet}/{replicate}/{location}_genome.fasta.all.gff'
#     output:
#         'from_cluster/maker_annotations/soap/{location}/{diet}/{replicate}/soap_{location}_{diet}_{replicate}.gff'
#     shell:
#         'mv {input} {output}'
#
# '''
# snakemake --cores 1 --use-conda \
# from_cluster/maker_annotations/soap/bristol/as/rep123/soap_bristol_as_rep123.gff
# '''

# rule rename_rep1_2_3:
#     input:
#         'from_cluster/maker_annotations/soap/{location}/{diet}/rep1_2_3/soap_{location}_{diet}_rep1_2_3.gff'
#     output:
#         'from_cluster/maker_annotations/soap/{location}/{diet}/rep123/soap_{location}_{diet}_rep123.gff'
#     shell:
#         'mv {input} {output}'
#
# '''
# snakemake --cores 1 \
# from_cluster/maker_annotations/soap/bristol/as/rep123/soap_bristol_as_rep123.gff
# '''


rule granges_liftover:
    input:
        script='scripts/annotations/granges.R',
        annotation='from_cluster/liftover_annotations/{location}/liftover_{location}.gtf'
    output:
        granges='output/annotations/granges/liftover/{location}/liftover_{location}.rds',
        granges_genes='output/annotations/granges_genes/liftover/{location}/liftover_{location}.rds'
    conda:
        '../envs/conda/bioconductor-rtracklayer=1.50.0.yaml'
    script:
        '../scripts/annotations/granges.R'

'''
snakemake --cores 1 --use-conda \
output/annotations/granges/liftover/altadena/liftover_altadena.rds
'''

rule granges_maker:
    ''' Make granges object from Maker annotation '''
    input:
        script='scripts/annotations/granges.R',
        annotation='from_cluster/maker_annotations/{annotation_type}/{location}/{diet}/{replicate}/{annotation_type}_{location}_{diet}_{replicate}.gff'
    output:
        granges='output/annotations/granges/{annotation_type}/{location}/{diet}/{replicate}/{annotation_type}_{location}_{diet}_{replicate}.rds',
        granges_genes='output/annotations/granges_genes/{annotation_type}/{location}/{diet}/{replicate}/{annotation_type}_{location}_{diet}_{replicate}.rds'
    conda:
        '../envs/conda/bioconductor-rtracklayer=1.50.0.yaml'
    script:
        '../scripts/annotations/granges.R'

'''
snakemake --cores 1 --use-conda -R \
output/annotations/granges/soap/altadena/all/rep123/soap_altadena_all_rep123.rds
'''

rule granges_stringtie:
    ''' Make granges object from Stringtie annotation '''
    input:
        script='scripts/annotations/granges.R',
        annotation='from_cluster/maker_annotations/{annotation_type}/{location}/{diet}/{replicate}/{annotation_type}_{location}_{diet}_{replicate}.gtf'
    output:
        granges='output/annotations/granges/{annotation_type}/{location}/{diet}/{replicate}/{annotation_type}_{location}_{diet}_{replicate}.rds',
        granges_genes='output/annotations/granges_genes/{annotation_type}/{location}/{diet}/{replicate}/{annotation_type}_{location}_{diet}_{replicate}.rds'
    conda:
        '../envs/conda/bioconductor-rtracklayer=1.50.0.yaml'
    script:
        '../scripts/annotations/granges.R'

'''
snakemake --cores 1 --use-conda \
output/annotations/granges/stringtie/bristol/as/rep2/stringtie_bristol_as_rep2.rds
'''

rule rename_reference:
    ''' Rename vc2010 reference annotation '''
    input:
        'input/wormbase_annotations/vc2010/c_elegans.PRJEB28388.WS279.annotations.gff3'
    output:
        'input/wormbase_annotations/vc2010/vc2010.gff3'
    shell:
        'mv {input} {output}'

rule granges_reference:
    ''' Make granges objects for reference annotation '''
    input:
        script='scripts/annotations/granges.R',
        annotation='input/wormbase_annotations/vc2010/vc2010.gff3'
    output:
        granges='output/annotations/granges/reference/vc2010/reference_vc2010.rds',
        granges_genes='output/annotations/granges_genes/reference/vc2010/reference_vc2010.rds'
    conda:
        '../envs/conda/bioconductor-rtracklayer=1.50.0.yaml'
    script:
        '../scripts/annotations/granges.R'

'''
snakemake --cores 1 --use-conda -R \
output/annotations/granges/liftover_bristol/liftover_bristol.rds
'''

rule coverage_liftover:
    ''' Calculatet coverage for liftover annotation '''
    input:
        script='scripts/annotations/coverage.R',
        granges='output/annotations/granges/{annotation_type}/{location}/{annotation_type}_{location}.rds'
    output:
        coverage='output/annotations/coverage/{annotation_type}/{location}/{annotation_type}_{location}.rds'
        #Also adds line to:
        #coverage_table='output/annotations/coverage/coverage_table.csv'
    conda:
        '../envs/conda/bioconductor-genomicranges=1.42.0.yaml'
    script:
        '../scripts/annotations/coverage.R'

'''
snakemake --cores 1 --use-conda -R \
output/annotations/coverage/liftover/bristol/liftover_bristol.rds
'''

rule coverage_new_annotations:
    ''' Calculate coverage for Maker or Stringtie annotation '''
    input:
        script='scripts/annotations/coverage.R',
        granges='output/annotations/granges/{annotation_type}/{location}/{diet}/{replicate}/{annotation_type}_{location}_{diet}_{replicate}.rds'
    output:
        coverage='output/annotations/coverage/{annotation_type}/{location}/{diet}/{replicate}/{annotation_type}_{location}_{diet}_{replicate}.rds'
        #Also adds line to:
        #coverage_table='output/annotations/coverage/coverage_table.txt'
    conda:
        '../envs/conda/bioconductor-genomicranges=1.42.0.yaml'
    script:
        '../scripts/annotations/coverage.R'

'''
snakemake --cores 1 --use-conda \
output/annotations/coverage/stringtie/altadena/as/rep123/stringtie_altadena_as_rep123.rds \
output/annotations/coverage/stringtie/altadena/bp/rep123/stringtie_altadena_bp_rep123.rds \
output/annotations/coverage/stringtie/altadena/hb101/rep123/stringtie_altadena_hb101_rep123.rds \
output/annotations/coverage/stringtie/altadena/m9/rep123/stringtie_altadena_m9_rep123.rds \
output/annotations/coverage/stringtie/altadena/op50/rep123/stringtie_altadena_op50_rep123.rds \
output/annotations/coverage/stringtie/altadena/pf/rep123/stringtie_altadena_pf_rep123.rds \
output/annotations/coverage/soap/altadena/as/rep123/soap_altadena_as_rep123.rds \
output/annotations/coverage/soap/altadena/bp/rep123/soap_altadena_bp_rep123.rds \
output/annotations/coverage/soap/altadena/hb101/rep123/soap_altadena_hb101_rep123.rds \
output/annotations/coverage/soap/altadena/m9/rep123/soap_altadena_m9_rep123.rds \
output/annotations/coverage/soap/altadena/op50/rep123/soap_altadena_op50_rep123.rds \
output/annotations/coverage/soap/altadena/pf/rep123/soap_altadena_pf_rep123.rds

snakemake --cores 1 --use-conda \
output/annotations/coverage/soap/bristol/allgff/rep123/soap_bristol_allgff_rep123.rds

'''

rule plot_coverage_comparison_bristol:
    ''' Compare Trinity, soap, stringtie, liftover and reference for Bristol '''
    input:
        script='scripts/annotations/plot_coverage_comparison_bristol.R',
        coverage_table='output/annotations/coverage/coverage_table.txt'
    output:
        coverage_plot='output/annotations/coverage/plots/coverage_plot_comparison_bristol.tiff'
    conda:
        '../envs/conda/r-ggplot2=3.3.1.yaml'
    script:
        '../scripts/annotations/plot_coverage_comparison_bristol.R'

'''
snakemake --cores 1 --use-conda -R \
output/annotations/coverage/plots/coverage_plot_comparison_bristol.tiff
'''

rule plot_coverage_comparison:
    ''' Compare soap, stringtie, liftover and reference for Bristol and Altadena '''
    input:
        script='scripts/annotations/plot_coverage_comparison.R',
        coverage_table='output/annotations/coverage/coverage_table.txt'
    output:
        coverage_plot='output/annotations/coverage/plots/coverage_plot_comparison.tiff'
    conda:
        '../envs/conda/r-ggplot2=3.3.1.yaml'
    script:
        '../scripts/annotations/plot_coverage_comparison.R'

'''
snakemake --cores 1 --use-conda output/annotations/coverage/plots/coverage_plot_comparison.tiff
'''

rule find_overlaps_maker:
    ''' Find and count overlaps between two Maker or Stringtie annotations '''
    input:
        script='scripts/annotations/find_overlaps.R',
        granges1='output/annotations/granges/{annotation_type1}/{location1}/{diet1}/{replicate1}/{annotation_type1}_{location1}_{diet1}_{replicate1}.rds',
        granges2='output/annotations/granges/{annotation_type2}/{location2}/{diet2}/{replicate2}/{annotation_type2}_{location2}_{diet2}_{replicate2}.rds'
    output:
        overlaps='output/annotations/overlaps/{annotation_type1}_{location1}_{diet1}_{replicate1}_vs_{annotation_type2}_{location2}_{diet2}_{replicate2}.rds',
        overlap_counts='output/annotations/overlap_counts/{annotation_type1}_{location1}_{diet1}_{replicate1}_vs_{annotation_type2}_{location2}_{diet2}_{replicate2}.rds'
    conda:
        '../envs/conda/bioconductor-genomicranges=1.42.0.yaml'
    script:
        '../scripts/annotations/find_overlaps.R'

'''
snakemake --cores 1 --use-conda \
output/annotations/overlaps/trinity_bristol_as_rep2_vs_soap_bristol_as_rep2.rds \
output/annotations/overlap_counts/trinity_bristol_as_rep2_vs_soap_bristol_as_rep2.rds
'''

rule find_overlaps_liftover_maker:
    ''' Find and count overlaps between liftover and new annotation '''
    input:
        script='scripts/annotations/find_overlaps.R',
        granges1='output/annotations/granges/liftover/{location1}/liftover_{location1}.rds',
        granges2='output/annotations/granges/{annotation_type2}/{location2}/{diet2}/{replicate2}/{annotation_type2}_{location2}_{diet2}_{replicate2}.rds'
    output:
        overlaps='output/annotations/overlaps/liftover_{location1}_vs_{annotation_type2}_{location2}_{diet2}_{replicate2}.rds',
        overlap_counts='output/annotations/overlap_counts/liftover_{location1}_vs_{annotation_type2}_{location2}_{diet2}_{replicate2}.rds'
    conda:
        '../envs/conda/bioconductor-genomicranges=1.42.0.yaml'
    script:
        '../scripts/annotations/find_overlaps.R'

'''
snakemake --cores 1 --use-conda \
output/annotations/overlaps/liftover_altadena_vs_soap_altadena_as_rep123.rds \
snakemake --cores 1 --use-conda \
output/annotations/overlaps/liftover_bristol_vs_soap_bristol_as_rep123.rds
'''

rule find_overlaps_liftover_maker_genes:
    ''' Fid overlaps between liftover and new annotation, genes only '''
    input:
        script='scripts/annotations/find_overlaps.R',
        granges1='output/annotations/granges_genes/liftover/{location1}/liftover_{location1}.rds',
        granges2='output/annotations/granges_genes/{annotation_type2}/{location2}/{diet2}/{replicate2}/{annotation_type2}_{location2}_{diet2}_{replicate2}.rds'
    output:
        overlaps='output/annotations/overlaps_genes/liftover_{location1}_vs_{annotation_type2}_{location2}_{diet2}_{replicate2}.rds',
        overlap_counts='output/annotations/overlap_counts_genes/liftover_{location1}_vs_{annotation_type2}_{location2}_{diet2}_{replicate2}.rds'
    conda:
        '../envs/conda/bioconductor-genomicranges=1.42.0.yaml'
    script:
        '../scripts/annotations/find_overlaps.R'

'''
snakemake --cores 1 --use-conda \
output/annotations/overlaps_genes/liftover_altadena_vs_soap_altadena_as_rep123.rds \
snakemake --cores 1 --use-conda \
output/annotations/overlaps/liftover_bristol_vs_soap_bristol_as_rep123.rds
'''

rule find_overlaps_liftover:
    ''' Find overlaps between two liftover annotations '''
    input:
        script='scripts/annotations/find_overlaps.R',
        granges1='output/annotations/granges/liftover/{location1}/liftover_{location1}.rds',
        granges2='output/annotations/granges/liftover/{location2}/liftover_{location2}.rds'
    output:
        overlaps='output/annotations/overlaps/liftover_{location1}_vs_liftover_{location2}.rds'
    conda:
        '../envs/conda/bioconductor-genomicranges=1.42.0.yaml'
    script:
        '../scripts/annotations/find_overlaps.R'

'''
snakemake --cores 1 --use-conda \
output/annotations/overlaps/liftover_bristol_vs_liftover_altadena.rds
'''

rule new_genes:
    ''' Make granges object containing genes found in new annotation but not liftover '''
    input:
        script='scripts/annotations/new_genes.R',
        overlap_counts_genes='output/annotations/overlap_counts_genes/liftover_{location}_vs_{annotation_type}_{location}_{diet}_{replicate}.rds',
        granges_genes='output/annotations/granges_genes/{annotation_type}/{location}/{diet}/{replicate}/{annotation_type}_{location}_{diet}_{replicate}.rds'
    output:
        granges_new_genes='output/annotations/granges_new_genes/{annotation_type}/{location}/{diet}/{replicate}/{annotation_type}_{location}_{diet}_{replicate}.rds'
        # Also adds line to output/annotations/new_genes/new_gene_counts.txt
    conda:
        '../envs/conda/bioconductor-plyranges=1.10.0.yaml'
    script:
        '../scripts/annotations/new_genes.R'

'''
snakemake --cores 1 --use-conda \
output/annotations/granges_new_genes/soap/bristol/as/rep123/soap_bristol_as_rep123.rds \
output/annotations/granges_new_genes/soap/bristol/bp/rep123/soap_bristol_bp_rep123.rds \
output/annotations/granges_new_genes/soap/bristol/hb101/rep123/soap_bristol_hb101_rep123.rds \
output/annotations/granges_new_genes/soap/bristol/m9/rep123/soap_bristol_m9_rep123.rds \
output/annotations/granges_new_genes/soap/bristol/op50/rep123/soap_bristol_op50_rep123.rds \
output/annotations/granges_new_genes/soap/bristol/pf/rep123/soap_bristol_pf_rep123.rds

snakemake --cores 1 --use-conda \
output/annotations/granges_new_genes/soap/bristol/all/rep123/soap_bristol_all_rep123.rds
'''

rule merge_new_genes:
    ''' Mak granges object containing only new genes found in all diet treatments '''
    input:
        script='scripts/annotations/merge_new_genes.R',
        granges_new_genes=expand('output/annotations/granges_new_genes/soap/{{location}}/{diet}/rep123/soap_{{location}}_{diet}_rep123.rds', \
        diet = DIETS)
    output:
        granges_new_genes='output/annotations/granges_new_genes_merged/{location}/{location}.rds'
    conda:
        '../envs/conda/bioconductor-genomicranges=1.42.0_bioconductor-plyranges=1.10.0.yaml'
    script:
        '../scripts/annotations/merge_new_genes.R'

'''
snakemake --cores 1 --use-conda \
output/annotations/granges_new_genes_merged/altadena/altadena.rds
'''

rule find_overlaps_liftover_all_maker_genes:
    ''' Finds and counts overlaps between all parts of liftover and just genes in new '''
    input:
        script='scripts/annotations/find_overlaps.R',
        granges1='output/annotations/granges/liftover/{location1}/liftover_{location1}.rds',
        granges2='output/annotations/granges_genes/{annotation_type2}/{location2}/{diet2}/{replicate2}/{annotation_type2}_{location2}_{diet2}_{replicate2}.rds'
    output:
        overlaps='output/annotations/overlaps_all_genes/liftover_{location1}_vs_{annotation_type2}_{location2}_{diet2}_{replicate2}.rds',
        overlap_counts='output/annotations/overlap_counts_all_genes/liftover_{location1}_vs_{annotation_type2}_{location2}_{diet2}_{replicate2}_genes.rds'
    conda:
        '../envs/conda/bioconductor-genomicranges=1.42.0.yaml'
    script:
        '../scripts/annotations/find_overlaps.R'

'''
snakemake --cores 1 --use-conda \
output/annotations/overlap_counts_all_genes/liftover_altadena_vs_soap_altadena_as_rep123_genes.rds
'''

rule new_genes2:
    input:
        script='scripts/annotations/new_genes.R',
        overlap_counts_genes='output/annotations/overlap_counts_all_genes/liftover_{location}_vs_{annotation_type}_{location}_{diet}_{replicate}_genes.rds',
        granges_genes='output/annotations/granges_genes/{annotation_type}/{location}/{diet}/{replicate}/{annotation_type}_{location}_{diet}_{replicate}.rds'
    output:
        granges_new_genes='output/annotations/granges_new_genes2/{annotation_type}/{location}/{diet}/{replicate}/{annotation_type}_{location}_{diet}_{replicate}.rds'
        # Also adds line to output/annotations/new_genes/new_gene_counts.txt
    conda:
        '../envs/conda/bioconductor-plyranges=1.10.0.yaml'
    script:
        '../scripts/annotations/new_genes.R'

'''
snakemake --cores 1 --use-conda \
output/annotations/granges_new_genes2/soap/bristol/as/rep123/soap_bristol_as_rep123.rds \
output/annotations/granges_new_genes/soap/bristol/bp/rep123/soap_bristol_bp_rep123.rds \
output/annotations/granges_new_genes/soap/bristol/hb101/rep123/soap_bristol_hb101_rep123.rds \
output/annotations/granges_new_genes/soap/bristol/m9/rep123/soap_bristol_m9_rep123.rds \
output/annotations/granges_new_genes/soap/bristol/op50/rep123/soap_bristol_op50_rep123.rds \
output/annotations/granges_new_genes/soap/bristol/pf/rep123/soap_bristol_pf_rep123.rds

snakemake --cores 1 --use-conda \
output/annotations/granges_new_genes/soap/bristol/all/rep123/soap_bristol_all_rep123.rds
'''

rule merge_new_genes2:
    input:
        script='scripts/annotations/merge_new_genes.R',
        granges_new_genes=expand('output/annotations/granges_new_genes2/soap/{{location}}/{diet}/rep123/soap_{{location}}_{diet}_rep123.rds', \
        diet = DIETS)
    output:
        granges_new_genes='output/annotations/granges_new_genes_merged2/{location}/{location}.rds'
    conda:
        '../envs/conda/bioconductor-genomicranges=1.42.0_bioconductor-plyranges=1.10.0.yaml'
    script:
        '../scripts/annotations/merge_new_genes.R'

'''
snakemake --cores 1 --use-conda \
output/annotations/granges_new_genes_merged2/altadena/altadena.rds
'''

rule new_genes3:
    input:
        script='scripts/annotations/new_genes.R',
        overlap_counts_genes='output/annotations/overlap_counts/liftover_{location}_vs_{annotation_type}_{location}_{diet}_{replicate}.rds',
        granges_genes='output/annotations/granges/{annotation_type}/{location}/{diet}/{replicate}/{annotation_type}_{location}_{diet}_{replicate}.rds'
    output:
        granges_new_genes='output/annotations/granges_new_genes3/{annotation_type}/{location}/{diet}/{replicate}/{annotation_type}_{location}_{diet}_{replicate}.rds'
        # Also adds line to output/annotations/new_genes/new_gene_counts.txt
    conda:
        '../envs/conda/bioconductor-plyranges=1.10.0.yaml'
    script:
        '../scripts/annotations/new_genes.R'

'''
snakemake --cores 1 --use-conda \
output/annotations/granges_new_genes3/soap/altadena/all/rep123/soap_altadena_all_rep123.rds
'''

rule find_sequences:
    '''Makes fasta files of all the new genes'''
    input:
        script='scripts/annotations/find_sequences.R',
        granges_new_genes='output/annotations/granges_new_genes_merged/{location}/{location}.rds',
        genome=expand('input/genomes/{{location}}_genome_chr/{{location}}_genome_{chromosome}.fa', \
        chromosome = CHROMOSOMES)
    output:
        sequences='output/annotations/new_gene_sequences/{location}/{location}.fa'
    conda:
        '../envs/conda/r-seqinr=3.4_5.yaml'
    script:
        '../scripts/annotations/find_sequences.R'

'''
snakemake --cores 1 --use-conda \
output/annotations/new_gene_sequences/altadena/altadena.fa
'''

rule blast_new_sequences_transcripts:
    ''' Blast new sequence against reference transcripsts '''
    input:
        sequences='output/annotations/new_gene_sequences/{location}/{location}.fa',
        reference_transcripts='input/reference_transcripts/c_elegans.PRJEB28388.WS279.mRNA_transcripts.fa'
    output:
        'output/annotations/blast_new_sequences/{location}/{location}_blast.txt'
    conda:
        '../envs/conda/blast=2.10.1.yaml'
    shell:
        'blastn -query {input.sequences} -subject {input.reference_transcripts} \
        -outfmt 7 -out {output}'

'''
snakemake --cores 1 --use-conda \
output/annotations/blast_new_sequences/bristol/bristol_blast.txt
'''

rule extract_no_hits:
    ''' Extract results with no blast hits from above '''
    input:
        'output/annotations/blast_new_sequences/{location}/{location}_blast.txt'
    output:
        'output/annotations/blast_new_sequences/{location}/{location}_blast_no_hits.txt'
    shell:
        "grep '\ 0\ hits\ found' -B 3 {input} > {output}"

'''
snakemake --cores 1 --use-conda \
output/annotations/blast_new_sequences/altadena/altadena_blast_no_hits.txt
'''

rule blast_new_sequences_genome:
    ''' Blast new sequence against reference genome '''
    input:
        sequences='output/annotations/new_gene_sequences/{location}/{location}.fa',
        reference_transcripts='input/genomes/vc2010_reference/c_elegans.PRJEB28388.WS279.genomic.fa'
    output:
        'output/annotations/blast_new_sequences_genome/{location}/{location}_blast.txt'
    conda:
        '../envs/conda/blast=2.10.1.yaml'
    shell:
        'blastn -query {input.sequences} -subject {input.reference_transcripts} \
        -outfmt 7 -out {output}'

'''
snakemake --cores 1 --use-conda \
output/annotations/blast_new_sequences_genome/altadena/altadena_blast.txt
'''
