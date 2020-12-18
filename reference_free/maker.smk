# running
rule maker:
    input:
        trinity='polyA/reference_free/trinity/{location}/{diet}/{replicate}/trinity/Trinity.fasta',
        genome='polyA/reference_free/repeatmasker/{location}/{location}_genome.fasta.masked'
    output:
        'polyA/reference_free/maker/{location}/{diet}/{replicate}/{location}_masked.all.gff'
    conda:
        '../envs/conda/maker=2.31.10.yaml'
    threads: 24
    shell:
        # move to correct directory
        'cd ~/polyA/reference_free/maker/{wildcards.location}/{wildcards.diet}/{wildcards.replicate}; '
        # produce control files for maker
        'maker -CTL; '
        # make adjustments to maker_opts.ctl file
        "cat maker_opts.ctl | "
        # set genome sequence
        "sed 's|^genome=|genome=../../../../repeatmasker/{wildcards.location}/{wildcards.location}_genome.fasta.masked|' | "
        # set transcriptome fasta file
        "sed 's|^est=|est=../../../../trinity/{wildcards.location}/{wildcards.diet}/{wildcards.replicate}/trinity/Trinity.fasta|' | "
        # disables Repeatmasker because licensing agreement with RepBase broken
        "sed 's|^model_org=all|model_org=|' | "
        # infer gene predictions from transcriptome
        "sed 's|^est2genome=0|est2genome=1|' > maker_opts2.ctl;"
        'mv maker_opts2.ctl maker_opts.ctl;'
        'maker -cpus {threads} -q; '
        'gff3_merge -d {wildcards.location}_masked.maker.output/{wildcards.location}_masked_master_datastore_index.log'

snakemake --cores 24 --use-conda --snakefile rules/maker.smk \
polyA/reference_free/maker/bristol/as/rep3/bristol_masked.all.gff

snakemake --cluster-config snakemake_profile/slurm.json --use-conda \
--profile snakemake_profile --cores 24 --snakefile rules/maker.smk \
polyA/reference_free/maker/bristol/as/rep3/bristol_masked.all.gff
