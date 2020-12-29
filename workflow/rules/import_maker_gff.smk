# on local machine
rule import_maker_gff:
    ''' Moves maker outputs from cluster to my local machine'''
    input:
        'sb2226@172.25.11.131:/polyA/reference_free/maker/{location}/{diet}/{included_replicates}/{location}_genome.fasta.all.gff'
    output:
        '/Users/Sarah/OneDrive/Documents/Uni/III/Project/Scripts/workflow/output/maker_gffs/{location}/{diet}/{location}_{diet}_{included_replicates}.gff'
    shell:
        'cd /Users/Sarah/OneDrive/Documents/Uni/III/Project/Scripts/workflow/output/maker_gffs/;'
        'scp {input} . ;'
        'mv {wildcards.location}_genome.fasta.all.gff {wildcards.location}_{wildcards.diet}_{wildcards.included_replicates}.gff'

'''
snakemake --use-conda --cores 1 --snakefile workflow/rules/import_maker_gff.smk \
/Users/Sarah/OneDrive/Documents/Uni/III/Project/Scripts/workflow/output/maker_gffs/bristol/as/bristol_as_rep3_2.gff
'''
