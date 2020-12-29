rule maker_replicates:
    input:
        premade='polyA/reference_free/maker/{location}/{diet}/rep3/{location}_genome.fasta.all.gff',
        genome='polyA/reference_free/repeatmasker/{location}/{location}_genome.fasta.masked',
        trinity='polyA/reference_free/trinity/{location}/{diet}/rep2/trinity/Trinity.fasta',
        reference_annotation='input/annotations/liftover/{location}/{location}_annotation.gtf'
    output:
        'polyA/reference_free/maker/{location}/{diet}/rep3_2/{location}_genome.fasta.all.gff'
    conda:
        '../envs/conda/maker=2.31.10.yaml'
    threads: 24
    shell:
        # move to correct directory
        'cd ~/polyA/reference_free/maker/{wildcards.location}/{wildcards.diet}/rep3_2; '
        # produce control files for maker
        'maker -CTL; '
        # make adjustments to maker_opts.ctl file
        "cat maker_opts.ctl | "
        # previously produced Maker annotation to be improved
        "sed 's|^maker_gff=|maker_gff=../rep3/{wildcards.location}_genome.fasta.all.gff|' | "
        # use ESTs in maker_gff
        "sed 's|^est_pass=0|est_pass=1|' | "
        # keeps names from previous run
        "sed 's|^model_pass=0|model_pass=1|' | "
        "sed 's|^map_forward=0|map_forward=1|' | "

        # same as previous
        # set genome sequence
        "sed 's|^genome=|genome=../../../../repeatmasker/{wildcards.location}/{wildcards.location}_genome.fasta.masked|' | "
        # set transcriptome fasta file - rep2
        "sed 's|^est=|est=../../../../trinity/{wildcards.location}/{wildcards.diet}/rep2/trinity/Trinity.fasta|' | "
        # uses names from existing liftover annotation where possible
        "sed 's|^model_gff=|model_gff=../../../../../../input/annotations/liftover/{wildcards.location}/{wildcards.location}_annotation.gtf|' | "
        # disables Repeatmasker because licensing agreement with RepBase broken
        "sed 's|^model_org=all|model_org=|' | "
        # infer gene predictions from transcriptome
        "sed 's|^est2genome=0|est2genome=1|' > maker_opts2.ctl;"
        'mv maker_opts2.ctl maker_opts.ctl;'
        'maker -cpus {threads} -q;'
        'gff3_merge -d {wildcards.location}_genome.fasta.maker.output/{wildcards.location}_genome.fasta_master_datastore_index.log'

snakemake --cluster-config snakemake_profile/slurm.json --use-conda \
--profile snakemake_profile --cores 24 --snakefile rules/maker_replicates.smk \
polyA/reference_free/maker/bristol/as/rep3_2/bristol_genome.fasta.all.gff

snakemake --cores 24 --use-conda --snakefile rules/maker_replicates.smk \
polyA/reference_free/maker/bristol/as/rep3_2/bristol_genome.fasta.all.gff
