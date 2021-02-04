# cat copy_gene_association.WS279.wb.c_elegans | sed '/^!/d' > GO.txt

rule wbid2go_mapping:
    input:
        script='scripts/wbid2go_mapping.R',
        go_terms='from_cluster/gene_ontology/GO.txt'
    output:
        wbid2go='output/topgo/wbid2go.txt'
    conda:
        '../envs/conda/r-tidyverse=1.2.1.yaml'
    script:
        '../scripts/wbid2go_mapping.R'

'''
snakemake --cores 1 --use-conda output/topgo/wbid2go.txt
'''
