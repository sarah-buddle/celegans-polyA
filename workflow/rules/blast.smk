''' Perform homology analysis of new genes using BLAST '''

LOCATIONS = ['bristol']
DIETS = ['as','bp','hb101','m9','op50','pf']
REPLICATES = ['rep1','rep2','rep3']
CHROMOSOMES = ['chrI','chrII','chrIII','chrIV','chrV','chrX']
SPECIES = ['brenneri.PRJNA20035','briggsae.PRJNA10731','inopinata.PRJDB5687', \
'latens.PRJNA248912','nigoni.PRJNA384657','remanei.PRJNA53967', \
'sinica.PRJNA194557','tropicalis.PRJNA53597']

rule find_sequences:
    '''Makes fasta files of all the new genes'''
    input:
        script='scripts/blast/find_sequences.R',
        allnewgenes=ancient('output/annotations/combined_annotations/combinedall/altadena/allnewgenes_altadena.rds'),
        genome=expand('input/genomes/{{location}}_genome_chr/{{location}}_genome_{chromosome}.fa', \
        chromosome = CHROMOSOMES)
    output:
        sequences='output/blast/new_gene_sequences/{location}/{location}.fa'
    conda:
        '../envs/conda/r-seqinr=3.4_5.yaml'
    script:
        '../scripts/blast/find_sequences.R'

'''
snakemake --cores 1 --use-conda -R \
output/blast/new_gene_sequences/bristol/bristol.fa
'''

rule blast_new_sequences_genome:
    ''' Blast new sequence against VC2010 reference genome '''
    input:
        sequences='output/blast/new_gene_sequences/{location}/{location}.fa',
        reference_genome='input/genomes/vc2010_reference/c_elegans.PRJEB28388.WS279.genomic.fa'
    output:
        'output/blast/vc2010_genome_hits/{location}/vc2010_genome_hits_{location}.txt'
    conda:
        '../envs/conda/blast=2.10.1.yaml'
    shell:
        'blastn -query {input.sequences} -subject {input.reference_genome} \
        -outfmt 7 -out {output}'

'''
snakemake --cores 1 --use-conda -R \
output/blast/vc2010_genome_hits/altadena/vc2010_genome_hits_altadena.txt
'''

rule extract_no_hits_genome:
    ''' Extract results with no blast hits from above '''
    input:
        'output/blast/vc2010_genome_hits/{location}/vc2010_genome_hits_{location}.txt'
    output:
        'output/blast/vc2010_genome_no_hits/{location}/vc2010_genome_no_hits_{location}.txt'
    shell:
        "grep '\ 0\ hits\ found' -B 3 {input} > {output}"

'''
snakemake --cores 1 --use-conda -R \
output/blast/vc2010_genome_no_hits/altadena/vc2010_genome_no_hits_altadena.txt
'''

rule no_hits_list_genome:
    ''' Extract gene names from above list '''
    input:
        'output/blast/vc2010_genome_no_hits/{location}/vc2010_genome_no_hits_{location}.txt'
    output:
        'output/blast/vc2010_genome_no_hits/{location}/names_vc2010_genome_no_hits_{location}.txt'
    shell:
        "grep 'Query' {input} | cut -c 10- > {output}"

'''
snakemake --cores 1 --use-conda \
output/blast/vc2010_genome_no_hits/altadena/names_vc2010_genome_no_hits_altadena.txt
'''

rule blast_new_sequences_transcripts:
    ''' Blast new sequence against reference transcripts '''
    input:
        sequences='output/blast/new_gene_sequences/{location}/{location}.fa',
        reference_transcripts='input/reference_transcripts/c_elegans.PRJEB28388.WS279.mRNA_transcripts.fa'
    output:
        'output/blast/vc2010_transcripts_hits/{location}/vc2010_transcripts_hits_{location}.txt'
    conda:
        '../envs/conda/blast=2.10.1.yaml'
    shell:
        'blastn -query {input.sequences} -subject {input.reference_transcripts} \
        -outfmt 7 -out {output}'

'''
snakemake --cores 1 --use-conda \
output/blast/vc2010_transcripts_hits/altadena/vc2010_transcripts_hits_altadena.txt
'''

rule extract_no_hits_transcripts:
    ''' Extract results with no blast hits from above '''
    input:
        'output/blast/vc2010_transcripts_hits/{location}/vc2010_transcripts_hits_{location}.txt'
    output:
        'output/blast/vc2010_transcripts_no_hits/{location}/vc2010_transcripts_no_hits_{location}.txt'
    shell:
        "grep '\ 0\ hits\ found' -B 3 {input} > {output}"

'''
snakemake --cores 1 --use-conda \
output/blast/vc2010_transcripts_no_hits/altadena/vc2010_transcripts_no_hits_altadena.txt
'''


rule no_hits_list_transcripts:
    ''' Extract gene names from above list '''
    input:
        'output/blast/vc2010_transcripts_no_hits/{location}/vc2010_transcripts_no_hits_{location}.txt'
    output:
        'output/blast/vc2010_transcripts_no_hits/{location}/names_vc2010_transcripts_no_hits_{location}.txt'
    shell:
        "grep 'Query' {input} | cut -c 10- > {output}"

'''
snakemake --cores 1 --use-conda \
output/blast/vc2010_transcripts_no_hits/altadena/names_vc2010_transcripts_no_hits_altadena.txt
'''

rule no_hits_sequences_transcripts:
    ''' Find sequences of genes with no hits to VC2010 transcripts '''
    input:
        script='scripts/blast/no_hits_sequences.R',
        no_hits_names='output/blast/vc2010_transcripts_no_hits/{location}/names_vc2010_transcripts_no_hits_{location}.txt',
        sequences='output/blast/new_gene_sequences/{location}/{location}.fa'
    output:
        no_hits='output/blast/vc2010_transcripts_no_hits/{location}/vc2010_transcripts_no_hits_{location}.fa'
    conda:
        '../envs/conda/r-seqinr=3.4_5.yaml'
    script:
        '../scripts/blast/no_hits_sequences.R'

'''
snakemake --cores 1 --use-conda \
output/blast/vc2010_transcripts_no_hits/altadena/vc2010_transcripts_no_hits_altadena.fa
'''

rule blast_new_sequences_transcripts_related_species2:
    ''' Blast new sequence against transcripts from related species '''
    input:
        sequences='output/blast/vc2010_transcripts_no_hits/{location}/vc2010_transcripts_no_hits_{location}.fa',
        reference_transcripts='input/reference_transcripts/c_{species}.WS279.mRNA_transcripts.fa'
    output:
        'output/blast/related_species_transcripts/{location}/c.{species}/{location}_c.{species}_transcripts.txt'
    conda:
        '../envs/conda/blast=2.10.1.yaml'
    shell:
        'blastn -query {input.sequences} -subject {input.reference_transcripts} \
        -outfmt 7 -out {output}'

'''
snakemake --cores 1 --use-conda \
output/blast/related_species_transcripts/altadena/c.brenneri.PRJNA20035/altadena_c.brenneri.PRJNA20035_transcripts.txt \
output/blast/related_species_transcripts/altadena/c.briggsae.PRJNA10731/altadena_c.briggsae.PRJNA10731_transcripts.txt \
output/blast/related_species_transcripts/altadena/c.inopinata.PRJDB5687/altadena_c.inopinata.PRJDB5687_transcripts.txt \
output/blast/related_species_transcripts/altadena/c.latens.PRJNA248912/altadena_c.latens.PRJNA248912_transcripts.txt \
output/blast/related_species_transcripts/altadena/c.nigoni.PRJNA384657/altadena_c.nigoni.PRJNA384657_transcripts.txt \
output/blast/related_species_transcripts/altadena/c.remanei.PRJNA53967/altadena_c.remanei.PRJNA53967_transcripts.txt \
output/blast/related_species_transcripts/altadena/c.sinica.PRJNA194557/altadena_c.sinica.PRJNA194557_transcripts.txt \
output/blast/related_species_transcripts/altadena/c.tropicalis.PRJNA53597/altadena_c.tropicalis.PRJNA53597_transcripts.txt
'''

rule extract_hits_related_species_transcripts:
    ''' Extract hits from above '''
    input:
        expand('output/blast/related_species_transcripts/{{location}}/c.{species}/{{location}}_c.{species}_transcripts.txt', \
        species = SPECIES)
    output:
        'output/blast/related_species_transcripts/{location}/related_species_transcripts_hits_{location}.txt'
    shell:
        "grep -r -v '#' {input} > {output}"

'''
snakemake --cores 1 --use-conda \
output/blast/related_species_transcripts/bristol/related_species_transcripts_hits_bristol.txt
'''

rule hits_related_species_list_transcripts:
    ''' Extract gene names from above '''
    input:
        'output/blast/related_species_transcripts/{location}/related_species_transcripts_hits_{location}.txt'
    output:
        'output/blast/related_species_transcripts/{location}/names_related_species_transcripts_hits_{location}.txt'
    shell:
        """cut -f1 {input} | cut -d ":" -f2 > {output}"""

'''
snakemake --cores 1 --use-conda \
output/blast/related_species_transcripts/altadena/names_related_species_transcripts_hits_altadena.txt
'''



rule blast_new_sequences_proteins:
    ''' Blast new sequence against reference proteins '''
    input:
        sequences='output/blast/new_gene_sequences/{location}/{location}.fa',
        reference_proteins='input/reference_proteins/c_elegans.PRJEB28388.WS279.protein.fa'
    output:
        'output/blast/vc2010_proteins_hits/{location}/vc2010_proteins_hits_{location}.txt'
    conda:
        '../envs/conda/blast=2.10.1.yaml'
    shell:
        'blastx -query {input.sequences} -subject {input.reference_proteins} \
        -outfmt 7 -out {output}'

'''
snakemake --cores 1 --use-conda \
output/blast/vc2010_proteins_hits/altadena/vc2010_proteins_hits_altadena.txt
'''

rule extract_no_hits_proteins:
    ''' Extract results with no blast hits from above '''
    input:
        'output/blast/vc2010_proteins_hits/{location}/vc2010_proteins_hits_{location}.txt'
    output:
        'output/blast/vc2010_proteins_no_hits/{location}/vc2010_proteins_no_hits_{location}.txt'
    shell:
        "grep '\ 0\ hits\ found' -B 3 {input} > {output}"

'''
snakemake --cores 1 --use-conda \
output/blast/vc2010_proteins_no_hits/altadena/vc2010_proteins_no_hits_altadena.txt
'''

rule vc2010_proteins_no_hits_list:
    ''' Extract gene names from above list '''
    input:
        'output/blast/vc2010_proteins_no_hits/{location}/vc2010_proteins_no_hits_{location}.txt'
    output:
        'output/blast/vc2010_proteins_no_hits/{location}/names_vc2010_proteins_no_hits_{location}.txt'
    shell:
        "grep 'Query' {input} | cut -c 10- > {output}"

'''
snakemake --cores 1 --use-conda \
output/blast/vc2010_proteins_no_hits/altadena/names_vc2010_proteins_no_hits_altadena.txt
'''

rule vc2010_proteins_no_hits_sequences:
    input:
        script='scripts/blast/no_hits_sequences.R',
        no_hits_names='output/blast/vc2010_proteins_no_hits/{location}/names_vc2010_proteins_no_hits_{location}.txt',
        sequences='output/blast/new_gene_sequences/{location}/{location}.fa'
    output:
        no_hits='output/blast/vc2010_proteins_no_hits/{location}/vc2010_proteins_no_hits_{location}.fa'
    conda:
        '../envs/conda/r-seqinr=3.4_5.yaml'
    script:
        '../scripts/blast/no_hits_sequences.R'

'''
snakemake --cores 1 --use-conda \
output/blast/vc2010_proteins_no_hits/altadena/vc2010_proteins_no_hits_altadena.fa
'''

rule blast_new_sequences_proteins_related_species2:
    ''' Blast new sequence against reference proteins '''
    input:
        sequences='output/blast/vc2010_proteins_no_hits/{location}/vc2010_proteins_no_hits_{location}.fa',
        reference_proteins='input/reference_proteins/c_{species}.WS279.protein.fa'
    output:
        'output/blast/related_species_proteins/{location}/c.{species}/{location}_c.{species}_proteins.txt'
    conda:
        '../envs/conda/blast=2.10.1.yaml'
    shell:
        'blastx -query {input.sequences} -subject {input.reference_proteins} \
        -outfmt 7 -out {output}'

'''
snakemake --cores 1 --use-conda \
output/blast/related_species_proteins/altadena/c.brenneri.PRJNA20035/altadena_c.brenneri.PRJNA20035_proteins.txt \
output/blast/related_species_proteins/altadena/c.briggsae.PRJNA10731/altadena_c.briggsae.PRJNA10731_proteins.txt \
output/blast/related_species_proteins/altadena/c.inopinata.PRJDB5687/altadena_c.inopinata.PRJDB5687_proteins.txt \
output/blast/related_species_proteins/altadena/c.latens.PRJNA248912/altadena_c.latens.PRJNA248912_proteins.txt \
output/blast/related_species_proteins/altadena/c.nigoni.PRJNA384657/altadena_c.nigoni.PRJNA384657_proteins.txt \
output/blast/related_species_proteins/altadena/c.remanei.PRJNA53967/altadena_c.remanei.PRJNA53967_proteins.txt \
output/blast/related_species_proteins/altadena/c.sinica.PRJNA194557/altadena_c.sinica.PRJNA194557_proteins.txt \
output/blast/related_species_proteins/altadena/c.tropicalis.PRJNA53597/altadena_c.tropicalis.PRJNA53597_proteins.txt
'''

rule extract_hits_related_species_proteins:
    ''' Extract hits from above '''
    input:
        expand('output/blast/related_species_proteins/{{location}}/c.{species}/{{location}}_c.{species}_proteins.txt', \
        species = SPECIES)
    output:
        'output/blast/related_species_proteins/{location}/related_species_proteins_hits_{location}.txt'
    shell:
        "grep -r -v '#' {input} > {output}"

'''
snakemake --cores 1 --use-conda \
output/blast/related_species_proteins/bristol/related_species_proteins_hits_bristol.txt
'''

rule hits_related_species_list_proteins:
    ''' Extract gene names from above '''
    input:
        'output/blast/related_species_proteins/{location}/related_species_proteins_hits_{location}.txt'
    output:
        'output/blast/related_species_proteins/{location}/names_related_species_proteins_hits_{location}.txt'
    shell:
        """cut -f1 {input} | cut -d ":" -f2 > {output}"""

'''
snakemake --cores 1 --use-conda \
output/blast/related_species_proteins/bristol/names_related_species_proteins_hits_bristol.txt
'''

# Bug in snakemake rule but works fine running script in R
rule venn:
    ''' Venn diagram to visualise above '''
    input:
        script='scripts/blast/blast_table.R',
        new_gene_counts='output/annotations/new_gene_expression/data/{location}/new_gene_counts_{location}.rds',
        no_hits_genome='output/blast/vc2010_genome_no_hits/{location}/names_vc2010_genome_no_hits_{location}.txt',
        no_hits_transcripts='output/blast/vc2010_transcripts_no_hits/{location}/names_vc2010_transcripts_no_hits_{location}.txt',
        no_hits_proteins='output/blast/vc2010_proteins_no_hits/{location}/names_vc2010_proteins_no_hits_{location}.txt',
        hits_related_species='output/blast/related_species_transcripts/{location}/names_related_species_transcripts_hits_{location}.txt'
    output:
        'output/blast/venn/{location}/venn_{location}.tiff'
    conda:
        '../envs/conda/r-tidyverse=1.2.1.yaml'
    script:
        '../scripts/blast/blast_table.R'

'''
snakemake --cores 1 --use-conda -R \
output/blast/venn/bristol/venn_bristol.tiff
'''

rule tree:
    ''' Produce phylogenetic tree of Caenorhabditis (uses Cristian's data) '''
    input:
        script='scripts/blast/tree.R',
        tree='input/tree/SpeciesTree_rooted_node_labels.txt'
    output:
        plot='output/blast/tree/tree_plot.tiff'
    conda:
        '../envs/conda/bioconductor-ggtree=2.4.1.yaml'
    script:
        '../scripts/blast/tree.R'

'''
snakemake --cores 1 --use-conda -R \
output/blast/tree/tree_plot.tiff
'''

rule related_species_plot:
    ''' Plot of genes that have homology to related species '''
    input:
        script='scripts/blast/related_species_plot.R',
        counts_bristol='output/blast/count_table/bristol/count_table_bristol.rds',
        counts_altadena='output/blast/count_table/altadena/count_table_altadena.rds'
    output:
        related_species_plot='output/blast/related_species_plot/related_species_plot.tiff'
    conda:
        '../envs/conda/r-tidyverse=1.2.1.yaml'
    script:
        '../scripts/blast/related_species_plot.R'

'''
snakemake --cores 1 --use-conda -R -n \
output/blast/related_species_plot/related_species_plot.tiff
'''
