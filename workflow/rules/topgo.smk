''' Gene ontology analysis of differentially expressed genes '''

''' Run deseq2.smk '''

'''
# Rename file containing GO term associations from WormBase
cat copy_gene_association.WS279.wb.c_elegans | sed '/^!/d' > GO.txt
'''

rule wbid2go_mapping:
    ''' Create file mapping WB IDs to GO terms '''
    input:
        script='scripts/topgo/wbid2go_mapping.R',
        go_terms='from_cluster/gene_ontology/GO.txt'
    output:
        wbid2go='output/topgo/data/wbid2go.txt'
    conda:
        '../envs/conda/r-tidyverse=1.2.1.yaml'
    script:
        '../scripts/topgo/wbid2go_mapping.R'

'''
snakemake --cores 1 --use-conda output/topgo/data/wbid2go.txt
'''

rule topgo_data:
    ''' Gene ontology analysis of differentially expressed genes '''
    input:
        script='scripts/topgo/topgo_data.R',
        wbid2go='output/topgo/data/wbid2go.txt',
        dds_res='output/deseq2/results/{location}/{location}_{diet1}_{diet2}.rds'
    output:
        topgo_data='output/topgo/data/{location}/{location}_{diet1}_{diet2}.rds'
    conda:
        '../envs/conda/bioconductor-topgo=2.42.0.yaml'
    script:
        '../scripts/topgo/topgo_data.R'

'''
snakemake --cores 1 --use-conda output/topgo_data/bristol/bristol_as_op50.rds
'''

rule go_enrichment_test:
    ''' Enrichment test in gene ontology analysis of differentially expressed genes '''
    input:
        script='scripts/topgo/go_enrichment_test.R',
        topgo_data='output/topgo/data/{location}/{location}_{diet1}_{diet2}.rds'
    output:
        topgo_result='output/topgo/results/{algorithm}/{statistic}/{location}/{location}_{diet1}_{diet2}.rds'
    conda:
        '../envs/conda/bioconductor-topgo=2.42.0.yaml'
    script:
        '../scripts/topgo/go_enrichment_test.R'

'''
snakemake --cores 1 --use-conda output/topgo/results/classic/fisher/bristol/bristol_as_op50.rds
'''

rule topgo_gentable:
    ''' Generate table of enrichment test results '''
    input:
        script='scripts/topgo/topgo_gentable.R',
        topgo_data='output/topgo/data/{location}/{location}_{diet1}_{diet2}.rds',
        topgo_result_classic_fisher='output/topgo/results/classic/fisher/{location}/{location}_{diet1}_{diet2}.rds',
        topgo_result_classic_ks='output/topgo/results/classic/ks/{location}/{location}_{diet1}_{diet2}.rds',
        topgo_result_elim_ks='output/topgo/results/elim/ks/{location}/{location}_{diet1}_{diet2}.rds'
    output:
        topgo_gentable='output/topgo/gentable/{orderBy}/{location}/{location}_{diet1}_{diet2}.rds'
    conda:
        '../envs/conda/bioconductor-topgo=2.42.0.yaml'
    script:
        '../scripts/topgo/topgo_gentable.R'

'''
snakemake --cores 1 --use-conda -R \
output/topgo/gentable/classicFisher/bristol/bristol_m9_op50.rds \
output/topgo/gentable/classicFisher/bristol/bristol_m9_pf.rds \
output/topgo/gentable/classicFisher/bristol/bristol_pf_op50.rds \
output/topgo/gentable/classicFisher/bristol/bristol_pf_as.rds \
output/topgo/gentable/classicFisher/altadena/altadena_m9_op50.rds \
output/topgo/gentable/classicFisher/altadena/altadena_m9_pf.rds \
output/topgo/gentable/classicFisher/altadena/altadena_pf_op50.rds \
output/topgo/gentable/classicFisher/altadena/altadena_pf_as.rds
'''

rule topgo_master_gentable1:
    ''' Comine gentables from multiple comparisons '''
    input:
        script='scripts/topgo/topgo_master_gentable.R',
        topgo_gentable=['output/topgo/gentable/classicFisher/bristol/bristol_m9_op50.rds', \
                'output/topgo/gentable/classicFisher/bristol/bristol_m9_pf.rds', \
                'output/topgo/gentable/classicFisher/bristol/bristol_pf_op50.rds', \
                'output/topgo/gentable/classicFisher/bristol/bristol_pf_as.rds', \
                'output/topgo/gentable/classicFisher/bristol/bristol_as_op50.rds', \
                'output/topgo/gentable/classicFisher/altadena/altadena_m9_op50.rds', \
                'output/topgo/gentable/classicFisher/altadena/altadena_m9_pf.rds', \
                'output/topgo/gentable/classicFisher/altadena/altadena_pf_op50.rds', \
                'output/topgo/gentable/classicFisher/altadena/altadena_pf_as.rds', \
                'output/topgo/gentable/classicFisher/altadena/altadena_as_op50.rds']
    output:
        master_gentable='output/topgo/master_gentable/1/master_gentable.rds'
    conda:
        '../envs/conda/r-tidyverse=1.2.1.yaml'
    script:
        '../scripts/topgo/topgo_master_gentable.R'

'''
snakemake --cores 1 --use-conda -R \
output/topgo/master_gentable/1/master_gentable.rds
'''

rule rrvgo:
    ''' Cluster GO terms by similarity '''
    input:
        script='scripts/topgo/rrvgo.R',
        master_gentable='output/topgo/master_gentable/1/master_gentable.rds'
    output:
        rrvgo='output/topgo/rvvgo/rrvgo.rds'
    conda:
        '../envs/conda/bioconductor-rrvgo=1.2.0_bioconductor-org.ce.eg.db=3.12.0.yaml'
    script:
        '../scripts/topgo/rrvgo.R'

rule topgo_plot:
    ''' Plot of GO analysis '''
    input:
        script='scripts/topgo/topgo_plot.R',
        master_gentable='output/topgo/master_gentable/1/master_gentable.rds',
        rrvgo='output/topgo/rvvgo/rrvgo.rds'
    output:
        topgo_plot='output/topgo/plots/1/topgo_plot_1.tiff'
    conda:
        '../envs/conda/r-tidyverse=1.2.1.yaml'
    script:
        '../scripts/topgo/topgo_plot.R'

'''
snakemake --cores 1 --use-conda -R --allowed-rules topgo_plot \
output/topgo/plots/1/topgo_plot_1.tiff
'''
