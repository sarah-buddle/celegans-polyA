rule maker_noest:
    input:
        genome='polyA/reference_free/repeatmasker/{location}/{location}_genome.fasta.masked',
        #trinity='polyA/reference_free/trinity/{location}/{diet}/{replicate}/trinity/Trinity.fasta',
        reference_annotation='input/annotations/liftover/{location}/{location}_annotation.gtf'
    output:
        'polyA/reference_free/maker/{location}/noest/{location}_genome.fasta.all.gff'
    conda:
        '../envs/conda/maker=2.31.10.yaml'
    threads: 8
    shell:
        # move to correct directory
        #'cd ~/polyA/reference_free/maker/{wildcards.location}/{wildcards.diet}/{wildcards.replicate}; '
        'cd polyA/reference_free/maker/{wildcards.location}/noest; '
        # produce control files for maker
        'maker -CTL; '
        # make adjustments to maker_opts.ctl file
        "cat maker_opts.ctl | "
        # set genome sequence
        "sed 's|^genome=|genome=../../../repeatmasker/{wildcards.location}/{wildcards.location}_genome.fasta.masked|' | "
        # set transcriptome fasta file
        #"sed 's|^est=|est=../../../../trinity/{wildcards.location}/{wildcards.diet}/{wildcards.replicate}/trinity/Trinity.fasta|' | "
        # uses names from existing liftover annotation where possible
        "sed 's|^model_gff=|model_gff=../../../../../input/annotations/liftover/{wildcards.location}/{wildcards.location}_annotation.gtf|' | "
        # disables Repeatmasker because licensing agreement with RepBase broken
        "sed 's|^model_org=all|model_org=|' > maker_opts2.ctl;"
        # infer gene predictions from transcriptome
        #"sed 's|^est2genome=0|est2genome=1|' > maker_opts2.ctl;"
        'mv maker_opts2.ctl maker_opts.ctl;'
        'maker -cpus {threads} -q;'
        'gff3_merge -d {wildcards.location}_genome.fasta.maker.output/{wildcards.location}_genome.fasta_master_datastore_index.log'

snakemake --cluster-config snakemake_profile/slurm.json --use-conda \
--profile snakemake_profile --cores 8 --snakefile rules/maker_noest.smk \
polyA/reference_free/maker/bristol/noest/bristol_genome.fasta.all.gff
