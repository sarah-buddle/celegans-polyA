LOCATIONS = ['bristol']
DIETS = ['as','bp','hb101','m9','op50','pf']
REPLICATES = ['rep1','rep2','rep3']

'''
scp -r sb2226@172.25.11.131://mnt/home1/miska/sb2226/workflow/output/polyA/reference_free/maker_soap_export \
from_cluster/maker_annotations/soap
'''

rule rename_maker:
    input:
        'from_cluster/maker_annotations/soap/maker_soap_export/{location}/{diet}/{replicate}/{location}_genome.fasta.all.gff'
    output:
        'from_cluster/maker_annotations/soap/{location}/{diet}/{replicate}/soap_{location}_{diet}_{replicate}.gff'
    shell:
        'mv {input} {output}'

'''
snakemake --cores 1 --use-conda \
from_cluster/maker_annotations/soap/bristol/as/rep1_2_3/soap_bristol_as_rep1_2_3.gff
'''

rule granges_liftover:
    input:
        script='scripts/annotations/granges.R',
        annotation='from_cluster/liftover_annotations/{location}/liftover_{location}.gtf'
    output:
        granges='output/annotations/granges/liftover/{location}/liftover_{location}.RData',
        granges_genes='output/annotations/granges_genes/liftover/{location}/liftover_{location}.RData'
    conda:
        '../envs/conda/bioconductor-rtracklayer=1.50.0.yaml'
    script:
        '../scripts/annotations/granges.R'

'''
snakemake --cores 1 --use-conda \
output/annotations/granges/liftover/altadena/liftover_altadena.RData
'''

rule granges_maker:
    input:
        script='scripts/annotations/granges.R',
        annotation='from_cluster/maker_annotations/{annotation_type}/{location}/{diet}/{replicate}/{annotation_type}_{location}_{diet}_{replicate}.gff'
    output:
        granges='output/annotations/granges/{annotation_type}/{location}/{diet}/{replicate}/{annotation_type}_{location}_{diet}_{replicate}.RData',
        granges_genes='output/annotations/granges_genes/{annotation_type}/{location}/{diet}/{replicate}/{annotation_type}_{location}_{diet}_{replicate}.RData'
    conda:
        '../envs/conda/bioconductor-rtracklayer=1.50.0.yaml'
    script:
        '../scripts/annotations/granges.R'

'''
snakemake --cores 1 --use-conda \
output/annotations/granges/trinity/bristol/as/rep2/trinity_bristol_as_rep2.RData
'''

rule coverage_liftover:
    input:
        script='scripts/annotations/coverage.R',
        granges='output/annotations/granges/{annotation_type}/{location}/{annotation_type}_{location}.RData'
    output:
        coverage='output/annotations/coverage/{annotation_type}/{location}/{annotation_type}_{location}.RData'
        #Also adds line to:
        #coverage_table='output/annotations/coverage/coverage_table.csv'
    conda:
        '../envs/conda/bioconductor-genomicranges=1.42.0.yaml'
    script:
        '../scripts/annotations/coverage.R'

'''
snakemake --cores 1 --use-conda \
output/annotations/coverage/liftover/bristol/liftover_bristol.RData
'''

rule coverage_maker:
    input:
        script='scripts/annotations/coverage.R',
        granges='output/annotations/granges/{annotation_type}/{location}/{diet}/{replicate}/{annotation_type}_{location}_{diet}_{replicate}.RData'
    output:
        coverage='output/annotations/coverage/{annotation_type}/{location}/{diet}/{replicate}/{annotation_type}_{location}_{diet}_{replicate}.RData'
        #Also adds line to:
        #coverage_table='output/annotations/coverage/coverage_table.txt'
    conda:
        '../envs/conda/bioconductor-genomicranges=1.42.0.yaml'
    script:
        '../scripts/annotations/coverage.R'

'''
snakemake --cores 1 --use-conda \
output/annotations/coverage/soap/bristol/as/rep1_2_3/soap_bristol_as_rep1_2_3.RData \
output/annotations/coverage/soap/bristol/bp/rep1_2_3/soap_bristol_bp_rep1_2_3.RData \
output/annotations/coverage/soap/bristol/hb101/rep1_2_3/soap_bristol_hb101_rep1_2_3.RData \
output/annotations/coverage/soap/bristol/m9/rep1_2_3/soap_bristol_m9_rep1_2_3.RData \
output/annotations/coverage/soap/bristol/op50/rep1_2_3/soap_bristol_op50_rep1_2_3.RData \
output/annotations/coverage/soap/bristol/pf/rep1_2_3/soap_bristol_pf_rep1_2_3.RData
'''

rule plot_coverage:
    input:
        script='scripts/annotations/plot_coverage.R',
        coverage_table='output/annotations/coverage/coverage_table.txt'
    output:
        coverage_plot='output/annotations/coverage/plots/coverage_plot.png'
    conda:
        '../envs/conda/r-ggplot2=3.3.1.yaml'
    script:
        '../scripts/annotations/plot_coverage.R'

'''
snakemake --cores 1 --use-conda output/annotations/coverage/plots/coverage_plot.png
'''

rule find_overlaps_maker:
    input:
        script='scripts/annotations/find_overlaps.R',
        granges1='output/annotations/granges/{annotation_type1}/{location1}/{diet1}/{replicate1}/{annotation_type1}_{location1}_{diet1}_{replicate1}.RData',
        granges2='output/annotations/granges/{annotation_type2}/{location2}/{diet2}/{replicate2}/{annotation_type2}_{location2}_{diet2}_{replicate2}.RData'
    output:
        overlaps='output/annotations/overlaps/{annotation_type1}_{location1}_{diet1}_{replicate1}_vs_{annotation_type2}_{location2}_{diet2}_{replicate2}.RData',
        overlap_counts='output/annotations/overlap_counts/{annotation_type1}_{location1}_{diet1}_{replicate1}_vs_{annotation_type2}_{location2}_{diet2}_{replicate2}.RData'
    conda:
        '../envs/conda/bioconductor-genomicranges=1.42.0.yaml'
    script:
        '../scripts/annotations/find_overlaps.R'

'''
snakemake --cores 1 --use-conda \
output/annotations/overlaps/trinity_bristol_as_rep2_vs_soap_bristol_as_rep2.RData \
output/annotations/overlap_counts/trinity_bristol_as_rep2_vs_soap_bristol_as_rep2.RData
'''

rule find_overlaps_liftover_maker:
    input:
        script='scripts/annotations/find_overlaps.R',
        granges1='output/annotations/granges/liftover/{location1}/liftover_{location1}.RData',
        granges2='output/annotations/granges/{annotation_type2}/{location2}/{diet2}/{replicate2}/{annotation_type2}_{location2}_{diet2}_{replicate2}.RData'
    output:
        overlaps='output/annotations/overlaps/liftover_{location1}_vs_{annotation_type2}_{location2}_{diet2}_{replicate2}.RData'
    conda:
        '../envs/conda/bioconductor-genomicranges=1.42.0.yaml'
    script:
        '../scripts/annotations/find_overlaps.R'

'''
snakemake --cores 1 --use-conda \
output/annotations/overlaps/liftover_bristol_vs_soap_bristol_as_rep2.RData
'''

rule find_overlaps_liftover:
    input:
        script='scripts/annotations/find_overlaps.R',
        granges1='output/annotations/granges/liftover/{location1}/liftover_{location1}.RData',
        granges2='output/annotations/granges/liftover/{location2}/liftover_{location2}.RData'
    output:
        overlaps='output/annotations/overlaps/liftover_{location1}_vs_liftover_{location2}.RData'
    conda:
        '../envs/conda/bioconductor-genomicranges=1.42.0.yaml'
    script:
        '../scripts/annotations/find_overlaps.R'

'''
snakemake --cores 1 --use-conda \
output/annotations/overlaps/liftover_bristol_vs_liftover_altadena.RData
'''
