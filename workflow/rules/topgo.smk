# cat copy_gene_association.WS279.wb.c_elegans | sed '/^!/d' > GO.txt

rule wbid2go_mapping:
    ''' Create file mapping WB IDs to GO terms '''
    input:
        script='scripts/topgo/wbid2go_mapping.R',
        go_terms='from_cluster/gene_ontology/GO.txt'
    output:
        wbid2go='output/topgo_data/wbid2go.txt'
    conda:
        '../envs/conda/r-tidyverse=1.2.1.yaml'
    script:
        '../scripts/topgo/wbid2go_mapping.R'

'''
snakemake --cores 1 --use-conda output/topgo/wbid2go.txt
'''

rule topgo_data:
    ''' Gene ontology analysis of differentially expressed genes '''
    input:
        script='scripts/topgo/topgo_data.R',
        wbid2go='output/topgo_data/wbid2go.txt',
        dds_res='output/deseq2_data/{location}/{location}_{diet1}_{diet2}.RData'
    output:
        topgo_data='output/topgo_data/{location}/{location}_{diet1}_{diet2}.RData'
    conda:
        '../envs/conda/bioconductor-topgo=2.42.0.yaml'
    script:
        '../scripts/topgo/topgo_data.R'

'''
snakemake --cores 1 --use-conda output/topgo_data/bristol/bristol_as_op50.RData
'''

rule go_enrichment_test:
    ''' Enrichment test in gene ontology analysis of differentially expressed genes '''
    input:
        script='scripts/topgo/go_enrichment_test.R',
        topgo_data='output/topgo_data/{location}/{location}_{diet1}_{diet2}.RData'
    output:
        topgo_result='output/topgo_result/{algorithm}/{statistic}/{location}/{location}_{diet1}_{diet2}.RData'
    conda:
        '../envs/conda/bioconductor-topgo=2.42.0.yaml'
    script:
        '../scripts/topgo/go_enrichment_test.R'

'''
snakemake --cores 1 --use-conda output/topgo_result/classic/fisher/bristol/bristol_as_op50.RData
'''

rule topgo_gentable:
    input:
        script='scripts/topgo/topgo_gentable.R',
        topgo_data='output/topgo_data/{location}/{location}_{diet1}_{diet2}.RData',
        topgo_result_classic_fisher='output/topgo_result/classic/fisher/{location}/{location}_{diet1}_{diet2}.RData',
        topgo_result_classic_ks='output/topgo_result/classic/ks/{location}/{location}_{diet1}_{diet2}.RData',
        topgo_result_elim_ks='output/topgo_result/elim/ks/{location}/{location}_{diet1}_{diet2}.RData'
    output:
        topgo_gentable='output/topgo_gentable/{orderBy}/{location}/{location}_{diet1}_{diet2}.RData'
    conda:
        '../envs/conda/bioconductor-topgo=2.42.0.yaml'
    script:
        '../scripts/topgo/topgo_gentable.R'

'''
snakemake --cores 1 --use-conda \
output/topgo_gentable/classicFisher/bristol/bristol_m9_op50.RData \
output/topgo_gentable/classicFisher/altadena/altadena_m9_op50.RData
'''
