''' De novo assemble transcripts using Trinity and run MAKER2 '''

''' NOt used in report due to slow speed '''

rule trinity:
    ''' De novo transcriptome assembly '''
    input:
        trimmed1="output/polyA/QC/trimmed_fastq/{location}/{diet}/{replicate}/reads1/{location}_{diet}_{replicate}_1_trimmed.fastq",
        trimmed2="output/polyA/QC/trimmed_fastq/{location}/{diet}/{replicate}/reads2/{location}_{diet}_{replicate}_2_trimmed.fastq"
    output:
        "output/polyA/reference_free/trinity/{location}/{diet}/{replicate}/trinity/Trinity.fasta"
    conda:
        "../envs/conda/Trinity=2.11.0.yaml"
    threads: 24
    shell:
        "Trinity --seqType fq --max_memory 60G --CPU {threads} \
        --left {input.trimmed1} --right {input.trimmed2} \
        --output output/polyA/reference_free/trinity/{wildcards.location}/{wildcards.diet}/{wildcards.replicate}/trinity/"

rule maker:
    ''' Run maker on single biological replicate '''
    input:
        genome='output/polyA/reference_free/repeatmasker/{location}/{location}_genome.fasta.masked',
        trinity='output/polyA/reference_free/trinity/{location}/{diet}/{replicate}/trinity/Trinity.fasta',
        reference_annotation='input/annotations/liftover/{location}/{location}_annotation.gtf'
    output:
        'output/polyA/reference_free/maker/{location}/{diet}/{replicate}/{location}_genome.fasta.all.gff'
    conda:
        '../envs/conda/maker=2.31.10.yaml'
    threads: 8
    shell:
        # move to correct directory
        'cd output/polyA/reference_free/maker/{wildcards.location}/{wildcards.diet}/{wildcards.replicate}; '
        # produce control files for maker
        'maker -CTL; '
        # make adjustments to maker_opts.ctl file
        "cat maker_opts.ctl | "
        # set genome sequence
        "sed 's|^genome=|genome=../../../../repeatmasker/{wildcards.location}/{wildcards.location}_genome.fasta.masked|' | "
        # set transcriptome fasta file
        "sed 's|^est=|est=../../../../trinity/{wildcards.location}/{wildcards.diet}/{wildcards.replicate}/trinity/Trinity.fasta|' | "
        # uses names from existing liftover annotation where possible
        "sed 's|^model_gff=|model_gff=../../../../../../../input/annotations/liftover/{wildcards.location}/{wildcards.location}_annotation.gtf|' | "
        # disables Repeatmasker because licensing agreement with RepBase broken
        "sed 's|^model_org=all|model_org=|' | "
        # infer gene predictions from transcriptome
        "sed 's|^est2genome=0|est2genome=1|' > maker_opts2.ctl;"
        'mv maker_opts2.ctl maker_opts.ctl;'
        'maker -cpus {threads} -q;'
        'gff3_merge -d {wildcards.location}_genome.fasta.maker.output/{wildcards.location}_genome.fasta_master_datastore_index.log'

'''
snakemake --profile ../snakemake_profile \
output/polyA/reference_free/maker/bristol/m9/rep2/bristol_genome.fasta.all.gff
'''


rule maker_rep2:
    ''' Run Maker on rep2 using info from rep1 '''
    input:
        premade='polyA/reference_free/maker/{location}/{diet}/rep1/{location}_genome.fasta.all.gff',
        genome='polyA/reference_free/repeatmasker/{location}/{location}_genome.fasta.masked',
        trinity='polyA/reference_free/trinity/{location}/{diet}/rep2/trinity/Trinity.fasta',
        reference_annotation='input/annotations/liftover/{location}/{location}_annotation.gtf'
    output:
        'polyA/reference_free/maker/{location}/{diet}/rep1_2/{location}_genome.fasta.all.gff'
    conda:
        '../envs/conda/maker=2.31.10.yaml'
    threads: 8
    shell:
        # move to correct directory
        'cd output/polyA/reference_free/maker/{wildcards.location}/{wildcards.diet}/rep1_2; '
        # produce control files for maker
        'maker -CTL; '
        # make adjustments to maker_opts.ctl file
        "cat maker_opts.ctl | "
        # previously produced Maker annotation to be improved
        "sed 's|^maker_gff=|maker_gff=../rep1/{wildcards.location}_genome.fasta.all.gff|' | "
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
        "sed 's|^model_gff=|model_gff=../../../../../../../input/annotations/liftover/{wildcards.location}/{wildcards.location}_annotation.gtf|' | "
        # disables Repeatmasker because licensing agreement with RepBase broken
        "sed 's|^model_org=all|model_org=|' | "
        # infer gene predictions from transcriptome
        "sed 's|^est2genome=0|est2genome=1|' > maker_opts2.ctl;"
        'mv maker_opts2.ctl maker_opts.ctl;'
        'maker -cpus {threads} -q;'
        'gff3_merge -d {wildcards.location}_genome.fasta.maker.output/{wildcards.location}_genome.fasta_master_datastore_index.log'

rule maker_rep3:
    ''' Run Maker on rep3 using info from rep1 and 2 '''
    input:
        premade='polyA/reference_free/maker/{location}/{diet}/rep1_2/{location}_genome.fasta.all.gff',
        genome='polyA/reference_free/repeatmasker/{location}/{location}_genome.fasta.masked',
        trinity='polyA/reference_free/trinity/{location}/{diet}/rep3/trinity/Trinity.fasta',
        reference_annotation='input/annotations/liftover/{location}/{location}_annotation.gtf'
    output:
        'polyA/reference_free/maker/{location}/{diet}/rep1_2_3/{location}_genome.fasta.all.gff'
    conda:
        '../envs/conda/maker=2.31.10.yaml'
    threads: 8
    shell:
        # move to correct directory
        'cd output/polyA/reference_free/maker/{wildcards.location}/{wildcards.diet}/rep1_2_3; '
        # produce control files for maker
        'maker -CTL; '
        # make adjustments to maker_opts.ctl file
        "cat maker_opts.ctl | "
        # previously produced Maker annotation to be improved
        "sed 's|^maker_gff=|maker_gff=../rep1_2/{wildcards.location}_genome.fasta.all.gff|' | "
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
        "sed 's|^model_gff=|model_gff=../../../../../../../input/annotations/liftover/{wildcards.location}/{wildcards.location}_annotation.gtf|' | "
        # disables Repeatmasker because licensing agreement with RepBase broken
        "sed 's|^model_org=all|model_org=|' | "
        # infer gene predictions from transcriptome
        "sed 's|^est2genome=0|est2genome=1|' > maker_opts2.ctl;"
        'mv maker_opts2.ctl maker_opts.ctl;'
        'maker -cpus {threads} -q;'
        'gff3_merge -d {wildcards.location}_genome.fasta.maker.output/{wildcards.location}_genome.fasta_master_datastore_index.log'
