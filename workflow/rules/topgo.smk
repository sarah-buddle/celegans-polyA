# cat copy_gene_association.WS279.wb.c_elegans | sed '/^!/d' > GO.txt

rule wbid2go_mapping:
    ''' Create file mapping WB IDs to GO terms '''
    input:
        script='scripts/wbid2go_mapping.R',
        go_terms='from_cluster/gene_ontology/GO.txt'
    output:
        wbid2go='output/topgo_data/wbid2go.txt'
    conda:
        '../envs/conda/r-tidyverse=1.2.1.yaml'
    script:
        '../scripts/wbid2go_mapping.R'

'''
snakemake --cores 1 --use-conda output/topgo/wbid2go.txt
'''

rule topgo_data:
    ''' Gene ontology analysis of differentially expressed genes '''
    input:
        script='scripts/topgo_data.R',
        wbid2go='output/topgo_data/wbid2go.txt',
        dds_res='output/deseq2_data/{location}/{location}_{diet1}_{diet2}.RData'
    output:
        topgo_data='output/topgo_data/{location}/{location}_{diet1}_{diet2}.RData'
    conda:
        '../envs/conda/bioconductor-topgo=2.42.0.yaml'
    script:
        '../scripts/topgo_data.R'

'''
snakemake --cores 1 --use-conda output/topgo_data/bristol/bristol_as_op50.RData
'''

rule go_enrichment_test:
    ''' Enrichment test in gene ontology analysis of differentially expressed genes '''
    input:
        script='scripts/go_enrichment_test.R',
        topgo_data='output/topgo_data/{location}/{location}_{diet1}_{diet2}.RData'
    output:
        topgo_result='output/topgo_result/{algorithm}/{statistic}/{location}/{location}_{diet1}_{diet2}.RData'
    conda:
        '../envs/conda/bioconductor-topgo=2.42.0.yaml'
    script:
        '../scripts/go_enrichment_test.R'

'''
snakemake --cores 1 --use-conda output/topgo_result/classic/fisher/bristol/bristol_as_op50.RData
'''
