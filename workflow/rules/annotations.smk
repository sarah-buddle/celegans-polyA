LOCATIONS = ['bristol']
DIETS = ['as','bp','hb101','m9','op50','pf']
REPLICATES = ['rep1','rep2','rep3']

'''
Run genome_guided.smk and soap_maker_all.smk on cluster
Export Stringtie, Maker and liftover annotations to workflow/from_cluster
'''

rule rename_maker:
    ''' Rename Maker annotations '''
    input:
        'from_cluster/maker_annotations/soap/{location}/{diet}/{replicate}/{location}_genome.fasta.all.gff'
    output:
        'from_cluster/maker_annotations/soap/{location}/{diet}/{replicate}/soap_{location}_{diet}_{replicate}.gff'
    shell:
        'mv {input} {output}'

'''
snakemake --cores 1 --use-conda \
from_cluster/maker_annotations/soap/bristol/as/rep123/soap_bristol_as_rep123.gff
'''

rule rename_rep1_2_3:
    ''' Rename rep1_2_3 to rep123 - needed for snakemake pipeline to work properly '''
    input:
        'from_cluster/maker_annotations/soap/{location}/{diet}/rep1_2_3/soap_{location}_{diet}_rep1_2_3.gff'
    output:
        'from_cluster/maker_annotations/soap/{location}/{diet}/rep123/soap_{location}_{diet}_rep123.gff'
    shell:
        'mv {input} {output}'

'''
snakemake --cores 1 \
from_cluster/maker_annotations/soap/bristol/as/rep123/soap_bristol_as_rep123.gff
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
    ''' Make granges object for reference annotation '''
    input:
        script='scripts/annotations/granges_reference.R',
        annotation='input/wormbase_annotations/vc2010/vc2010.gff3'
    output:
        granges='output/annotations/granges/reference/vc2010/reference_vc2010.rds',
        granges_genes='output/annotations/granges_genes/reference/vc2010/reference_vc2010.rds'
    conda:
        '../envs/conda/bioconductor-rtracklayer=1.50.0.yaml'
    script:
        '../scripts/annotations/granges_reference.R'

'''
snakemake --cores 1 --use-conda -R \
output/annotations/granges_genes/reference/vc2010/reference_vc2010.rds
'''

rule granges_liftover:
    ''' Make granges objects for reference annotation '''
    input:
        script='scripts/annotations/granges_liftover.R',
        annotation='from_cluster/liftover_annotations/{location}/liftover_{location}.gtf'
    output:
        granges='output/annotations/granges/liftover/{location}/liftover_{location}.rds',
        granges_genes='output/annotations/granges_genes/liftover/{location}/liftover_{location}.rds'
    conda:
        '../envs/conda/bioconductor-rtracklayer=1.50.0.yaml'
    script:
        '../scripts/annotations/granges_liftover.R'

'''
snakemake --cores 1 --use-conda -R \
output/annotations/granges/liftover/altadena/liftover_altadena.rds \
output/annotations/granges/liftover/bristol/liftover_bristol.rds
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

rule combine_annotations_maker:
    ''' Combine Maker annotations from different diets and identify new genes'''
    input:
        script='scripts/annotations/combine_annotations_maker.R',
        granges_liftover='output/annotations/granges/liftover/{location}/liftover_{location}.rds',
        granges_soap=expand('output/annotations/granges_genes/soap/{{location}}/{diet}/rep123/soap_{{location}}_{diet}_rep123.rds', \
        diet = DIETS)
    output:
        combined_annotation_full='output/annotations/combined_annotations/combinedmaker/{location}/combined_annotation_full_{location}.rds',
        combined_annotation='output/annotations/combined_annotations/combinedmaker/{location}/combinedmakeronly_{location}.rds',
        combined_annotation_gtf='output/annotations/combined_annotations/combinedmakeronly/{location}/combinedmakeronly_{location}.gtf',
        new_genes='output/annotations/combined_annotations/combinedmaker/{location}/new_genes_{location}.rds',
        maker_combined_annotation='output/annotations/combined_annotations/combinedmaker/{location}/combinedmaker_{location}.gtf'
    conda:
        '../envs/conda/bioconductor-rtracklayer=1.50.0_bioconductor-genomicranges=1.42.0_bioconductor-plyranges=1.10.0.yaml'
    script:
        '../scripts/annotations/combine_annotations_maker.R'

'''
snakemake --cores 1 --use-conda \
output/annotations/combined_annotations/combinedmaker/bristol/combinedmaker_bristol.gtf
'''

rule combine_annotations_stringtie:
    ''' Identify new genes in StringtTie annotation '''
    input:
        script='scripts/annotations/combine_annotations_stringtie.R',
        granges_liftover='output/annotations/granges/liftover/{location}/liftover_{location}.rds',
        granges_stringtie='output/annotations/granges/stringtie/{location}/all/rep123/stringtie_{location}_all_rep123.rds'
    output:
        new_genes='output/annotations/combined_annotations/combinedstringtie/{location}/new_genes_{location}.rds',
        stringtie_combined_annotation='output/annotations/combined_annotations/combinedstringtie/{location}/combinedstringtie_{location}.gtf'
    conda:
        '../envs/conda/bioconductor-rtracklayer=1.50.0_bioconductor-genomicranges=1.42.0_bioconductor-plyranges=1.10.0.yaml'
    script:
        '../scripts/annotations/combine_annotations_stringtie.R'

'''
snakemake --cores 1 --use-conda \
output/annotations/combined_annotations/combinedstringtie/bristol/combinedstringtie_bristol.gtf
'''

rule combine_annotations_all:
    ''' Produce combined annotations for differential expression analysis and prodice list of new genes '''
    input:
        script='scripts/annotations/combine_annotations_all.R',
        liftover_genes='output/annotations/granges_genes/liftover/{location}/liftover_{location}.rds',
        maker_new_genes='output/annotations/combined_annotations/combinedmaker/{location}/new_genes_{location}.rds',
        stringtie_new_genes='output/annotations/combined_annotations/combinedstringtie/{location}/new_genes_{location}.rds'
    output:
        all_combined_annotation='output/annotations/combined_annotations/combinedall/{location}/combinedall_{location}.gtf',
        all_new_genes='output/annotations/combined_annotations/combinedall/{location}/allnewgenes_{location}.rds'
    conda:
        '../envs/conda/bioconductor-rtracklayer=1.50.0_bioconductor-genomicranges=1.42.0_bioconductor-plyranges=1.10.0.yaml'
    script:
        '../scripts/annotations/combine_annotations_all.R'

'''
snakemake --cores 1 --use-conda -R \
output/annotations/combined_annotations/combinedall/bristol/combinedall_bristol.gtf \
output/annotations/combined_annotations/combinedall/altadena/combinedall_altadena.gtf
'''

'''
Export the combined annotation to cluster workflow/input/annotations/combinedall
Run genome_guided_combined on cluster
'''

rule move_combined_annotation:
    ''' Move combined annotation to IGV directory for viewing '''
    input:
        'output/annotations/combined_annotations/{source}/{location}/{source}_{location}.gtf'
    output:
        'output/igv/gff/no_contig/{source}/{location}/all/rep123/{source}_{location}_all_rep123.gtf'
    shell:
        'cp {input} {output}'

'''
snakemake --cores 1 --use-conda \
output/igv/gff/no_contig/combinedall/bristol/all/rep123/combinedall_bristol_all_rep123.gtf
'''

rule coverage_reference_liftover:
    ''' Calculate coverage for reference or liftover annotation '''
    input:
        script='scripts/annotations/coverage.R',
        granges='output/annotations/granges_genes/{annotation_type}/{location}/{annotation_type}_{location}.rds'
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
output/annotations/coverage/reference/vc2010/reference_vc2010.rds \
output/annotations/coverage/liftover/bristol/liftover_bristol.rds \
output/annotations/coverage/liftover/altadena/liftover_altadena.rds
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
'''

rule coverage_combined_annotations:
    ''' Calculate coverage of combined annotations '''
    input:
        script='scripts/annotations/coverage_combined.R',
        annotation='output/annotations/combined_annotations/{source}/{location}/{source}_{location}.gtf'
    output:
        coverage='output/annotations/coverage/{source}/{location}/{source}_{location}.rds'
        #Also adds line to:
        #coverage_table='output/annotations/coverage/coverage_table.csv'
    conda:
        '../envs/conda/bioconductor-rtracklayer=1.50.0_bioconductor-genomicranges=1.42.0_bioconductor-plyranges=1.10.0.yaml'
    script:
        '../scripts/annotations/coverage_combined.R'

'''
snakemake --cores 1 --use-conda -R \
output/annotations/coverage/combinedmakeronly/altadena/combinedmakeronly_altadena.rds \
output/annotations/coverage/combinedmakeronly/bristol/combinedmakeronly_bristol.rds \
output/annotations/coverage/combinedmaker/altadena/combinedmaker_altadena.rds \
output/annotations/coverage/combinedstringtie/altadena/combinedstringtie_altadena.rds \
output/annotations/coverage/combinedall/altadena/combinedall_altadena.rds \
output/annotations/coverage/combinedmaker/bristol/combinedmaker_bristol.rds \
output/annotations/coverage/combinedstringtie/bristol/combinedstringtie_bristol.rds \
output/annotations/coverage/combinedall/bristol/combinedall_bristol.rds
'''

rule plot_coverage_comparison:
    ''' Produce plot comparing soap, stringtie, liftover and reference for Bristol and Altadena '''
    input:
        script='scripts/annotations/plot_coverage_comparison4.R',
        coverage_table='output/annotations/coverage/coverage_table.txt'
    output:
        coverage_plot_tiff='output/annotations/coverage/plots/coverage_plot_comparison.tiff',
        coverage_plot_rds='output/annotations/coverage/plots/coverage_plot_comparison.rds'
    conda:
        '../envs/conda/r-ggplot2=3.3.1.yaml'
    script:
        '../scripts/annotations/plot_coverage_comparison4.R'

'''
snakemake --cores 1 -R --use-conda \
output/annotations/coverage/plots/coverage_plot_comparison.tiff
'''

# Script isn't written for snakemake yet
rule compare_gene_counts:
    ''' Calculate and plot gene counts in different annotation types '''
    input:
        script='scripts/annotations/compare_gene_counts.R',
        granges_genes=['output/annotations/granges_genes/reference/vc2010/reference_vc2010.rds', \
        'output/annotations/granges_genes/liftover/altadena/liftover_altadena.rds', \
        'output/annotations/granges_genes/liftover/bristol/liftover_bristol.rds', \
        'output/annotations/granges/stringtie/altadena/all/rep123/stringtie_altadena_all_rep123.rds', \
        'output/annotations/granges/stringtie/bristol/all/rep123/stringtie_bristol_all_rep123.rds'],
        annotations=['output/annotations/combined_annotations/combinedall/altadena/combinedall_altadena.gtf', \
        'output/annotations/combined_annotations/combinedmakeronly/altadena/combinedmakeronly_altadena.gtf', \
        'output/annotations/combined_annotations/combinedall/bristol/combinedall_bristol.gtf', \
        'output/annotations/combined_annotations/combinedmakeronly/bristol/combinedmakeronly_bristol.gtf']
    output:
        gene_count_plot='output/annotations/gene_counts/gene_counts_plot.tiff',
        gene_count_plot_rds='output/annotations/gene_counts/gene_counts_plot.rds',
        gene_counts='output/annotations/gene_counts/gene_counts.tsv'
    conda:
        '../envs/conda/bioconductor-rtracklayer=1.50.0_r-ggplot2=3.3.1_bioconductor-plyranges=1.10.0.yaml'
    script:
        '../scripts/annotations/compare_gene_counts.R'

'''
snakemake --cores 1 -R --use-conda \
output/annotations/gene_counts/gene_counts_plot.tiff
'''

rule plot_coverage_genecounts:
    ''' Combine coverage and gene counts plots with shared legend '''
    input:
        coverage_plot='output/annotations/coverage/plots/coverage_plot_comparison.rds',
        gene_plot='output/annotations/gene_counts/gene_counts_plot.rds'
    output:
        coverage_genecounts_plot='output/annotations/coverage_genecounts_plot/coverage_genecounts_plot.tiff'
    conda:
        '../r-ggplot2=3.3.1_rggpubr=0.4.0_r-gridextra=2.3.yaml'
    script:
        '../scripts/annotations/plot_coverage_genecounts.R'

'''
snakemake --cores 1 --use-conda -R \
output/annotations/coverage_genecounts_plot/coverage_genecounts_plot.tiff
'''

'''
Export the combined annotation to cluster workflow/input/annotations/combinedall
Run genome_guided_combined on cluster
'''
