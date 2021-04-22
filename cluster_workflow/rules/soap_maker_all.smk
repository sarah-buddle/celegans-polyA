''' De novo assemble transcripts using SOAPdenovo-trans and run MAKER2 '''

''' When running soap, need to specify global file path so adjust the lines specified '''

rule soap:
    ''' De novo transcriptome assembly '''
    input:
        q1='output/polyA/QC/trimmed_fastq/{location}/{diet}/{replicate}/reads1/{location}_{diet}_{replicate}_1_trimmed.fastq',
        q2='output/polyA/QC/trimmed_fastq/{location}/{diet}/{replicate}/reads2/{location}_{diet}_{replicate}_2_trimmed.fastq'
    output:
        'output/polyA/reference_free/soap/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.scafSeq'
    conda:
        '../envs/conda/soapdenovo-trans=1.04.yaml'
    shell:
        'cd output/polyA/reference_free/soap/{wildcards.location}/{wildcards.diet}/{wildcards.replicate};'
        'cp ../../../../../../../scripts/soap_config.txt .;'
        "cat soap_config.txt | "
        "sed 's|^q1=|q1=/mnt/home1/miska/sb2226/workflow/{input.q1}|' | " # adjust global file path
        "sed 's|^q2=|q2=/mnt/home1/miska/sb2226/workflow/{input.q2}|' > soap_config2.txt;" # adjust global file path
        'mv soap_config2.txt soap_config.txt;'
        'SOAPdenovo-Trans-127mer all -s soap_config.txt -o {wildcards.location}_{wildcards.diet}_{wildcards.replicate}'

rule maker_soap:
    ''' Run maker on single biological replicate '''
    input:
        genome='output/polyA/reference_free/repeatmasker/{location}/{location}_genome.fasta.masked',
        soap='output/polyA/reference_free/soap/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.scafSeq',
        reference_annotation='input/annotations/liftover/{location}/{location}_annotation.gtf'
    output:
        'output/polyA/reference_free/maker_soap/{location}/{diet}/{replicate}/{location}_genome.fasta.all.gff'
    conda:
        '../envs/conda/maker=2.31.10.yaml'
    shell:
        # move to correct directory
        'cd output/polyA/reference_free/maker_soap/{wildcards.location}/{wildcards.diet}/{wildcards.replicate}; '
        # produce control files for maker
        'maker -CTL; '
        # make adjustments to maker_opts.ctl file
        "cat maker_opts.ctl | "
        # set genome sequence
        "sed 's|^genome=|genome=../../../../repeatmasker/{wildcards.location}/{wildcards.location}_genome.fasta.masked|' | "
        # set transcriptome fasta file
        "sed 's|^est=|est=../../../../soap/{wildcards.location}/{wildcards.diet}/{wildcards.replicate}/{wildcards.location}_{wildcards.diet}_{wildcards.replicate}.scafSeq|' | "
        # uses names from existing liftover annotation where possible
        "sed 's|^model_gff=|model_gff=../../../../../../../input/annotations/liftover/{wildcards.location}/{wildcards.location}_annotation.gtf|' | "
        # disables Repeatmasker because licensing agreement with RepBase broken
        "sed 's|^model_org=all|model_org=|' | "
        # infer gene predictions from transcriptome
        "sed 's|^est2genome=0|est2genome=1|' > maker_opts2.ctl;"
        'mv maker_opts2.ctl maker_opts.ctl;'
        'maker -q;'
        'gff3_merge -d {wildcards.location}_genome.fasta.maker.output/{wildcards.location}_genome.fasta_master_datastore_index.log'

'''
snakemake --profile ../snakemake_profile \
output/polyA/reference_free/maker_soap/bristol/as/rep2/bristol_genome.fasta.all.gff \
output/polyA/reference_free/maker_soap/bristol/bp/rep2/bristol_genome.fasta.all.gff \
output/polyA/reference_free/maker_soap/bristol/hb101/rep2/bristol_genome.fasta.all.gff \
output/polyA/reference_free/maker_soap/bristol/m9/rep2/bristol_genome.fasta.all.gff \
output/polyA/reference_free/maker_soap/bristol/op50/rep2/bristol_genome.fasta.all.gff \
output/polyA/reference_free/maker_soap/bristol/pf/rep2/bristol_genome.fasta.all.gff
output/polyA/reference_free/maker_soap/bristol/as/rep3/bristol_genome.fasta.all.gff \
output/polyA/reference_free/maker_soap/bristol/bp/rep3/bristol_genome.fasta.all.gff \
output/polyA/reference_free/maker_soap/bristol/hb101/rep3/bristol_genome.fasta.all.gff \
output/polyA/reference_free/maker_soap/bristol/m9/rep3/bristol_genome.fasta.all.gff \
output/polyA/reference_free/maker_soap/bristol/op50/rep3/bristol_genome.fasta.all.gff \
output/polyA/reference_free/maker_soap/bristol/pf/rep3/bristol_genome.fasta.all.gff
'''

rule maker_rep2_soap:
    ''' Run Maker on rep2 using info from rep1 '''
    input:
        premade='output/polyA/reference_free/maker_soap/{location}/{diet}/rep1/{location}_genome.fasta.all.gff',
        genome='output/polyA/reference_free/repeatmasker/{location}/{location}_genome.fasta.masked',
        soap='output/polyA/reference_free/soap/{location}/{diet}/rep2/{location}_{diet}_rep2.scafSeq',
        reference_annotation='input/annotations/liftover/{location}/{location}_annotation.gtf'
    output:
        'output/polyA/reference_free/maker_soap/{location}/{diet}/rep1_2/{location}_genome.fasta.all.gff'
    conda:
        '../envs/conda/maker=2.31.10.yaml'
    threads: 8
    shell:
        # move to correct directory
        'cd output/polyA/reference_free/maker_soap/{wildcards.location}/{wildcards.diet}/rep1_2; '
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
        "sed 's|^est=|est=../../../../soap/{wildcards.location}/{wildcards.diet}/rep2/{wildcards.location}_{wildcards.diet}_rep2.scafSeq|' | "
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
output/polyA/reference_free/maker_soap/altadena/as/rep1_2/altadena_genome.fasta.all.gff \
output/polyA/reference_free/maker_soap/altadena/bp/rep1_2/altadena_genome.fasta.all.gff \
output/polyA/reference_free/maker_soap/altadena/hb101/rep1_2/altadena_genome.fasta.all.gff \
output/polyA/reference_free/maker_soap/altadena/m9/rep1_2/altadena_genome.fasta.all.gff \
output/polyA/reference_free/maker_soap/altadena/op50/rep1_2/altadena_genome.fasta.all.gff \
output/polyA/reference_free/maker_soap/altadena/pf/rep1_2/altadena_genome.fasta.all.gff
'''

rule maker_rep3_soap:
    ''' Run Maker on rep3 using info from rep1 and 2 '''
    input:
        premade='output/polyA/reference_free/maker_soap/{location}/{diet}/rep1_2/{location}_genome.fasta.all.gff',
        genome='output/polyA/reference_free/repeatmasker/{location}/{location}_genome.fasta.masked',
        soap='output/polyA/reference_free/soap/{location}/{diet}/rep3/{location}_{diet}_rep3.scafSeq',
        reference_annotation='input/annotations/liftover/{location}/{location}_annotation.gtf'
    output:
        'output/polyA/reference_free/maker_soap/{location}/{diet}/rep1_2_3/{location}_genome.fasta.all.gff'
    conda:
        '../envs/conda/maker=2.31.10.yaml'
    threads: 8
    shell:
        # move to correct directory
        'cd output/polyA/reference_free/maker_soap/{wildcards.location}/{wildcards.diet}/rep1_2_3; '
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
        "sed 's|^est=|est=../../../../soap/{wildcards.location}/{wildcards.diet}/rep3/{wildcards.location}_{wildcards.diet}_rep3.scafSeq|' | "
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
output/polyA/reference_free/maker_soap/altadena/as/rep1_2_3/altadena_genome.fasta.all.gff \
output/polyA/reference_free/maker_soap/altadena/bp/rep1_2_3/altadena_genome.fasta.all.gff \
output/polyA/reference_free/maker_soap/altadena/hb101/rep1_2_3/altadena_genome.fasta.all.gff \
output/polyA/reference_free/maker_soap/altadena/m9/rep1_2_3/altadena_genome.fasta.all.gff \
output/polyA/reference_free/maker_soap/altadena/op50/rep1_2_3/altadena_genome.fasta.all.gff \
output/polyA/reference_free/maker_soap/altadena/pf/rep1_2_3/altadena_genome.fasta.all.gff
'''

rule maker_soap_export:
    ''' Move to single file for export '''
    input:
        'output/polyA/reference_free/maker_soap/{location}/{diet}/rep1_2_3/{location}_genome.fasta.all.gff'
    output:
        'output/polyA/reference_free/maker_soap_export/{location}/{diet}/rep1_2_3/{location}_genome.fasta.all.gff'
    shell:
        'cp {input} {output}'

'''
snakemake --profile ../snakemake_profile \
output/polyA/reference_free/maker_soap_export/altadena/as/rep1_2_3/altadena_genome.fasta.all.gff \
output/polyA/reference_free/maker_soap_export/altadena/bp/rep1_2_3/altadena_genome.fasta.all.gff \
output/polyA/reference_free/maker_soap_export/altadena/hb101/rep1_2_3/altadena_genome.fasta.all.gff \
output/polyA/reference_free/maker_soap_export/altadena/m9/rep1_2_3/altadena_genome.fasta.all.gff \
output/polyA/reference_free/maker_soap_export/altadena/op50/rep1_2_3/altadena_genome.fasta.all.gff \
output/polyA/reference_free/maker_soap_export/altadena/pf/rep1_2_3/altadena_genome.fasta.all.gff
'''

'''
Export output/polyA/reference_free/maker_soap_export to workflow/from_cluster on local machine
'''

''' Merge Maker-derived annotations - not used in final report because this step is done using a custom scrip in R '''

rule merge_maker_all_stringtie:
    ''' Merge maker annotations using stringtie-merge '''
    input:
        expand('output/polyA/reference_free/maker_soap/{{location}}/{diet}/rep1_2_3/{{location}}_genome.fasta.all.gff', \
        diet = DIETS)
    output:
        'output/polyA/reference_free/maker_soap_export/{location}/all/rep123/soap_{location}_all_rep123.gtf'
    conda:
        '../envs/conda/stringtie=2.1.4.yaml'
    shell:
        'stringtie --merge {input} -o {output}'

'''
snakemake --profile ../snakemake_profile \
output/polyA/reference_free/maker_soap_export/bristol/all/rep123/soap_bristol_all_rep123.gtf
'''

rule merge_maker_all_gffcompare:
    ''' Merge maker annotations using gffcompare '''
    input:
        expand('output/polyA/reference_free/maker_soap/{{location}}/{diet}/rep1_2_3/{{location}}_genome.fasta.all.gff', \
        diet = DIETS)
    output:
        'output/polyA/reference_free/maker_soap_export/{location}/allgff/rep123/soap_{location}_allgff_rep123.combined.gff'
    conda:
        '../envs/conda/gffcompare=0.11.2.yaml'
    shell:
        'gffcompare {input} -o output/polyA/reference_free/maker_soap_export/{wildcards.location}/allgff/rep123/soap_{wildcards.location}_allgff_rep123'

rule merge_maker_liftover_gffcompare:
    ''' Merge maker and liftover annotations using gffcompare '''
    input:
        maker=expand('output/polyA/reference_free/maker_soap/{{location}}/{diet}/rep1_2_3/{{location}}_genome.fasta.all.gff', diet = DIETS),
        liftover='input/annotations/liftover/{location}/{location}_annotation.gtf'
    output:
        'output/polyA/reference_free/maker_soap/{location}/allgffliftover/rep123/soap_{location}_allgffliftover_rep123.combined.gtf'
    conda:
        '../envs/conda/gffcompare=0.11.2.yaml'
    shell:
        'gffcompare {input.maker} -r {input.liftover} -o output/polyA/reference_free/maker_soap/{wildcards.location}/allgffliftover/rep123/soap_{wildcards.location}_allgffliftover_rep123'


'''
snakemake --profile ../snakemake_profile \
output/polyA/reference_free/maker_soap/bristol/allgffliftover/rep123/soap_bristol_allgffliftover_rep123.combined.gtf
'''
