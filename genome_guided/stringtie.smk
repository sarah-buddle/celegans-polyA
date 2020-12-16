rule stringtie:
    input:
        "polyA/alignment/sorted_bam/{sample}.bam"
    output:
        "polyA/assembly/stringtie/{sample}.gtf"
    conda:
        "../envs/conda/stringtie=2.1.4.yaml"
    shell:
        "stringtie {input} -o {output} -p 8 -G annotations/bristol_annotation.gtf -e"

snakemake --use-conda --cores 8 --snakefile rules/stringtie.smk \
polyA/assembly/stringtie/{as_rep1,as_rep2,as_rep3,\
bp_rep1,bp_rep2,bp_rep3,\
hb101_rep1,hb101_rep2,hb101_rep3,\
m9_rep1,m9_rep2,m9_rep3,\
op50_rep1,op50_rep2,op50_rep3,\
pf_rep1,pf_rep2,pf_rep3}.gtf
