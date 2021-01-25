rule soap:
    input:
        q1='output/polyA/QC/trimmed_fastq/{location}/{diet}/{replicate}/reads1/{location}_{diet}_{replicate}_1_trimmed.fastq',
        q2='output/polyA/QC/trimmed_fastq/{location}/{diet}/{replicate}/reads2/{location}_{diet}_{replicate}_2_trimmed.fastq'
    output:
        'output/polyA/reference_free/soap/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.scafSeq'
    conda:
        '../envs/conda/soapdenovo-trans=1.04.yaml'
    shell:
        'cd output/polyA/reference_free/soap/{wildcards.location}/{wildcards.diet}/{wildcards.replicate}/;'
        'cp ../../../../../../../scripts/soap_config.txt .;'
        "cat soap_config.txt | "
        "sed 's|^q1=|q1=/mnt/home1/miska/sb2226/workflow/{input.q1}|' | "
        "sed 's|^q2=|q2=/mnt/home1/miska/sb2226/workflow/{input.q2}|' > soap_config2.txt;"
        'mv soap_config2.txt soap_config.txt;'
        'SOAPdenovo-Trans-127mer all -s soap_config.txt -o {wildcards.location}_{wildcards.diet}_{wildcards.replicate}'

'''
snakemake --cluster-config ../snakemake_profile/slurm.json --use-conda \
--profile ../snakemake_profile --cores 8 --snakefile rules/soap.smk \
output/polyA/reference_free/soap/altadena/as/rep2/altadena_as_rep2.scafSeq

snakemake --cluster-config ../snakemake_profile/slurm.json --use-conda \
--profile ../snakemake_profile --cores 24 --snakefile rules/soap.smk \
output/polyA/reference_free/soap/altadena/as/rep3/altadena_as_rep3.scafSeq \
output/polyA/reference_free/soap/altadena/bp/rep1/altadena_bp_rep1.scafSeq \
output/polyA/reference_free/soap/altadena/bp/rep2/altadena_bp_rep2.scafSeq \
output/polyA/reference_free/soap/altadena/bp/rep3/altadena_bp_rep3.scafSeq \
output/polyA/reference_free/soap/altadena/hb101/rep1/altadena_hb101_rep1.scafSeq \
output/polyA/reference_free/soap/altadena/hb101/rep2/altadena_hb101_rep2.scafSeq \
output/polyA/reference_free/soap/altadena/hb101/rep3/altadena_hb101_rep3.scafSeq \
output/polyA/reference_free/soap/altadena/m9/rep1/altadena_m9_rep1.scafSeq \
output/polyA/reference_free/soap/altadena/m9/rep2/altadena_m9_rep2.scafSeq \
output/polyA/reference_free/soap/altadena/m9/rep3/altadena_m9_rep3.scafSeq \
output/polyA/reference_free/soap/altadena/op50/rep1/altadena_op50_rep1.scafSeq \
output/polyA/reference_free/soap/altadena/op50/rep2/altadena_op50_rep2.scafSeq \
output/polyA/reference_free/soap/altadena/op50/rep3/altadena_op50_rep3.scafSeq \
output/polyA/reference_free/soap/altadena/pf/rep1/altadena_pf_rep1.scafSeq \
output/polyA/reference_free/soap/altadena/pf/rep2/altadena_pf_rep2.scafSeq \
output/polyA/reference_free/soap/altadena/pf/rep3/altadena_pf_rep3.scafSeq
'''
