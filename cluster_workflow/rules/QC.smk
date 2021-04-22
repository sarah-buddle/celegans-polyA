LOCATIONS = ['altadena']
DIETS = ['as','bp','hb101','m9','op50','pf']
REPLICATES = ['rep1','rep2','rep3']
READS = ['1','2']

'''
# rename replicate directories to rep
find . -type d -name *replicate* | xargs rename 's/replicate/rep/' .
'''

rule rename:
    ''' Rename to follow conventions in remaining rules '''
    input:
        'input/polyA/{location}/{diet}/{replicate}/reads.cram'
    output:
        'output/polyA/QC/cram/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.cram'
    shell:
        'cp {input} {output}'

rule cramtobam:
    ''' Convert cram files to bam '''
    input:
        'output/polyA/QC/cram/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.cram'
    output:
        'output/polyA/QC/bam/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.bam'
    conda:
        '../envs/conda/samtools=1.11.yaml'
    shell:
        'samtools view -o {output} {input}'

rule bam_to_fastq:
    ''' Convert bam files to fastq '''
    input:
        'output/polyA/QC/bam/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.bam'
    output:
        'output/polyA/QC/unsorted_fastq/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}_unsorted.fastq'
    conda:
        '../envs/conda/samtools=1.11.yaml'
    shell:
        'samtools bam2fq {input} > {output}'

rule sort_fastq:
    ''' Sort fastq files '''
    input:
        'output/polyA/QC/unsorted_fastq/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}_unsorted.fastq'
    output:
        'output/polyA/QC/sorted_fastq/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}_sorted.fastq'
    conda:
        '../envs/conda/fastq-tools=0.8.3.yaml'
    shell:
        'fastq-sort {input} > {output}'

rule deinterleave:
    ''' Separate paired end reads '''
    input:
        'output/polyA/QC/sorted_fastq/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}_sorted.fastq'
    output:
        reads1 = 'output/polyA/QC/deinterleaved_fastq/{location}/{diet}/{replicate}/reads1/{location}_{diet}_{replicate}_1.fastq',
        reads2 = 'output/polyA/QC/deinterleaved_fastq/{location}/{diet}/{replicate}/reads2/{location}_{diet}_{replicate}_2.fastq'
    conda:
        '../envs/conda/seqtk=1.3.yaml'
    shell:
        'seqtk seq -1 {input} > {output.reads1}; seqtk seq -2 {input} > {output.reads2}'

rule fastqc_untrimmed:
    ''' Perform fastqc to check quality of untrimmed reads '''
    input:
        'output/polyA/QC/deinterleaved_fastq/{location}/{diet}/{replicate}/reads{reads}/{location}_{diet}_{replicate}_{reads}.fastq'
    output:
        html='output/polyA/QC/fastqc_untrimmed/reports/{location}/{diet}/{replicate}/reads{reads}/{location}_{diet}_{replicate}_{reads}_fastqc.html',
        zip='output/polyA/QC/fastqc_untrimmed/zip/{location}/{diet}/{replicate}/reads{reads}/{location}_{diet}_{replicate}_{reads}_fastqc.zip'
    conda:
        '../envs/conda/fastqc=0.11.9.yaml'
    shell:
        'fastqc {input} --outdir output/polyA/QC/fastqc_untrimmed/tempdir;'
        'mv output/polyA/QC/fastqc_untrimmed/tempdir/{wildcards.location}_{wildcards.diet}_{wildcards.replicate}_{wildcards.reads}_fastqc.html {output.html};'
        'mv output/polyA/QC/fastqc_untrimmed/tempdir/{wildcards.location}_{wildcards.diet}_{wildcards.replicate}_{wildcards.reads}_fastqc.zip {output.zip}'

rule multiqc_untrimmed:
    ''' Combine fastqc reports on untrimmed reads '''
    input:
        expand('output/polyA/QC/fastqc_untrimmed/zip/{location}/{diet}/{replicate}/reads{reads}/{location}_{diet}_{replicate}_{reads}_fastqc.zip', \
        location=LOCATIONS, diet=DIETS, replicate=REPLICATES, reads=READS)
    output:
        report='output/polyA/QC/multiqc_untrimmed/multiqc_untrimmed.html'
    conda:
        '../envs/conda/multiqc=1.9.yaml'
    shell:
        'multiqc {input} -o output/polyA/QC/multiqc_untrimmed/ -n multiqc_untrimmed.html'

rule cutadapt:
    ''' Trims adaptors and low quality sequences, discards reads shorter than 40nt'''
    input:
        reads1='output/polyA/QC/deinterleaved_fastq/{location}/{diet}/{replicate}/reads1/{location}_{diet}_{replicate}_1.fastq',
        reads2='output/polyA/QC/deinterleaved_fastq/{location}/{diet}/{replicate}/reads2/{location}_{diet}_{replicate}_2.fastq'
    output:
        trimmed1='output/polyA/QC/trimmed_fastq/{location}/{diet}/{replicate}/reads1/{location}_{diet}_{replicate}_1_trimmed.fastq',
        trimmed2='output/polyA/QC/trimmed_fastq/{location}/{diet}/{replicate}/reads2/{location}_{diet}_{replicate}_2_trimmed.fastq'
    #group: 'QC_group'
    shell:
        "cutadapt -q 20 --minimum-length 40 \
        -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACAAGACTGTATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAA \
        -A GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTTAGTCCACGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAA \
        -o {output.trimmed1} -p {output.trimmed2} \
        {input.reads1} {input.reads2}"

rule fastqc_trimmed:
    ''' Perform fastqc to check quality of trimmed reads '''
    input:
        'output/polyA/QC/trimmed_fastq/{location}/{diet}/{replicate}/reads{reads}/{location}_{diet}_{replicate}_{reads}_trimmed.fastq'
    output:
        html='output/polyA/QC/fastqc_trimmed/reports/{location}/{diet}/{replicate}/reads{reads}/{location}_{diet}_{replicate}_{reads}_trimmed_fastqc.html',
        zip='output/polyA/QC/fastqc_trimmed/zip/{location}/{diet}/{replicate}/reads{reads}/{location}_{diet}_{replicate}_{reads}_trimmed_fastqc.zip'
    conda:
        '../envs/conda/fastqc=0.11.9.yaml'
    shell:
        'fastqc {input} --outdir output/polyA/QC/fastqc_trimmed/;'
        'mv output/polyA/QC/fastqc_trimmed/{wildcards.location}_{wildcards.diet}_{wildcards.replicate}_{wildcards.reads}_trimmed_fastqc.html {output.html};'
        'mv output/polyA/QC/fastqc_trimmed/{wildcards.location}_{wildcards.diet}_{wildcards.replicate}_{wildcards.reads}_trimmed_fastqc.zip {output.zip}'

rule multiqc_trimmed:
    ''' Combine fastqc reports on trimmed reads '''
    input:
        expand('output/polyA/QC/fastqc_trimmed/zip/{location}/{diet}/{replicate}/reads{reads}/{location}_{diet}_{replicate}_{reads}_trimmed_fastqc.zip', \
        location=LOCATIONS, diet=DIETS, replicate=REPLICATES, reads=READS)
    output:
        report='output/polyA/QC/multiqc_trimmed/multiqc_trimmed.html'
    conda:
        '../envs/conda/multiqc=1.9.yaml'
    shell:
        'multiqc {input} -o output/polyA/QC/multiqc_trimmed/ -n multiqc_trimmed.html'


'''
snakemake --profile ../snakemake_profile \
output/polyA/QC/multiqc_untrimmed/multiqc_untrimmed.html \
output/polyA/QC/multiqc_trimmed/multiqc_trimmed.html
'''

'''
Export reports to local machine to view
'''
