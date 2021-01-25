rule maker_soap:
    input:
        genome='polyA/reference_free/repeatmasker/{location}/{location}_genome.fasta.masked',
        soap='polyA/reference_free/soap/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.scafSeq',
        reference_annotation='input/annotations/liftover/{location}/{location}_annotation.gtf'
    output:
        'polyA/reference_free/maker_soap/{location}/{diet}/{replicate}/{location}_genome.fasta.all.gff'
    conda:
        '../envs/conda/maker=2.31.10.yaml'
    shell:
        # move to correct directory
        'cd ~/polyA/reference_free/maker_soap/{wildcards.location}/{wildcards.diet}/{wildcards.replicate}; '
        # produce control files for maker
        'maker -CTL; '
        # make adjustments to maker_opts.ctl file
        "cat maker_opts.ctl | "
        # set genome sequence
        "sed 's|^genome=|genome=../../../../repeatmasker/{wildcards.location}/{wildcards.location}_genome.fasta.masked|' | "
        # set transcriptome fasta file
        "sed 's|^est=|est=../../../../soap/{wildcards.location}/{wildcards.diet}/{wildcards.replicate}/{wildcards.location}_{wildcards.diet}_{wildcards.replicate}.scafSeq|' | "
        # uses names from existing liftover annotation where possible
        "sed 's|^model_gff=|model_gff=../../../../../../input/annotations/liftover/{wildcards.location}/{wildcards.location}_annotation.gtf|' | "
        # disables Repeatmasker because licensing agreement with RepBase broken
        "sed 's|^model_org=all|model_org=|' | "
        # infer gene predictions from transcriptome
        "sed 's|^est2genome=0|est2genome=1|' > maker_opts2.ctl;"
        'mv maker_opts2.ctl maker_opts.ctl;'
        'maker -q;'
        'gff3_merge -d {wildcards.location}_genome.fasta.maker.output/{wildcards.location}_genome.fasta_master_datastore_index.log'

snakemake --cores 24 --use-conda --snakefile rules/maker.smk \
polyA/reference_free/maker/bristol/as/rep3/bristol_genome.fasta.all.gff

snakemake --cluster-config snakemake_profile/slurm.json --use-conda \
--profile snakemake_profile --cores 24 --snakefile rules/maker_soap.smk \
polyA/reference_free/maker_soap/bristol/{as,bp,hb101,m9,op50,pf}/rep1/bristol_genome.fasta.all.gff

scp sb2226@172.25.11.131:/mnt/home1/miska/sb2226/polyA/reference_free/maker_soap/bristol/as/rep2/bristol_genome.fasta.all.gff .
mv bristol_genome.fasta.all.gff maker_soap_bristol_as_rep2.gff

scp sb2226@172.25.11.131:/mnt/home1/miska/sb2226/polyA/reference_free/maker/bristol/as/rep2/bristol_genome.fasta.all.gff .
mv bristol_genome.fasta.all.gff maker_bristol_as_rep2.gff
