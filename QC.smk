location = ['bristol']
diet = ['as','bp','hb101','m9','op50','pf']
replicate = ['rep1','rep2','rep3']
reads = ['1','2']

rule bam_to_fastqc:
    ''' Convert bam files to fastq '''
    input:
        'polyA/QC/bam/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}.bam'
    output:
        'polyA/QC/unsorted_fastq/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}_unsorted.fastq'
    conda:
        '../envs/conda/samtools=1.11.yaml'
    shell:
        'samtools bam2fq {input} > {output}'

rule sort_fastq:
    ''' Sort fastq files '''
    input:
        'polyA/QC/unsorted_fastq/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}_unsorted.fastq'
    output:
        'polyA/QC/sorted_fastq/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}_sorted.fastq'
    conda:
        '../envs/conda/fastq-tools=0.8.3.yaml'
    shell:
        'sort-fastq {input} > {output}'

rule deinterleave:
    ''' Separate paired end reads '''
    input:
        'polyA/QC/sorted_fastq/{location}/{diet}/{replicate}/{location}_{diet}_{replicate}_sorted.fastq'
    output:
        reads1 = 'polyA/QC/deinterleaved_fastq/{location}/{diet}/{replicate}/reads1/{location}_{diet}_{replicate}_1.fastq',
        reads2 = 'polyA/QC/deinterleaved_fastq/{location}/{diet}/{replicate}/reads2/{location}_{diet}_{replicate}_2.fastq'
    conda:
        '../envs/conda/seqtk=1.3.yaml'
    shell:
        'seqtk seq -1 {input} > {output.reads1}; seqtk seq -2 {input} > {output.reads2}'

rule fastqc_untrimmed:
    ''' Perform fastqc to check quality of untrimmed reads '''
    input:
        'polyA/QC/deinterleaved_fastq/{location}/{diet}/{replicate}/reads{reads}/{location}_{diet}_{replicate}_{reads}.fastq'
    output:
        html='polyA/QC/fastqc_untrimmed/{location}/{diet}/{replicate}/reads{reads}/{location}_{diet}_{replicate}_{reads}.html'
        zip='polyA/QC/fastqc_untrimmed/{location}/{diet}/{replicate}/reads{reads}/{location}_{diet}_{replicate}_{reads}.zip'
    conda:
        '../envs/conda/fastqc=0.11.9.yaml'
    shell:
        'fastqc {input} -o polyA/QC/fastqc_untrimmed'

rule multiqc_untrimmed:
    ''' Combine fastqc reports on untrimmed reads '''
    input:
        'polyA/QC/fastqc_untrimmed'
    output:
        data='polyA/QC/multiqc_untrimmed/untrimmed_multiqc_data'
        report='polyA/QC/multiqc_untrimmed/untrimmed_multiqc.html'
    conda:
        '../envs/conda/multiqc=1.9.yaml'
    shell:
        'multiqc {input}'

rule export_QC_reports:
    ''' Put QC reports into a single folder for export to local machine for viewing '''
    input:
        fastqc_html='polyA/QC/fastqc_untrimmed/{location}/{diet}/{replicate}/reads{reads}/{location}_{diet}_{replicate}_{reads}.html'
        fastqc_zip='polyA/QC/fastqc_untrimmed/{location}/{diet}/{replicate}/reads{reads}/{location}_{diet}_{replicate}_{reads}.zip'

rule cutadapt:
    ''' Trims adaptors and low quality sequences, discards reads shorter than 40nt'''
    input:
        reads1='polyA/QC/deinterleaved_fastq/{location}/{diet}/{replicate}/reads1/{location}_{diet}_{replicate}_1.fastq'
        reads2='polyA/QC/deinterleaved_fastq/{location}/{diet}/{replicate}/reads2/{location}_{diet}_{replicate}_2.fastq'
    output:
        trimmed1='polyA/QC/trimmed_fastq/{location}/{diet}/{replicate}/reads1/{location}_{diet}_{replicate}_1_trimmed.fastq'
        trimmed2='polyA/QC/trimmed_fastq/{location}/{diet}/{replicate}/reads2/{location}_{diet}_{replicate}_2_trimmed.fastq'
    shell:
        "cutadapt -q 20 --minimum-length 40 \
        -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACAAGACTGTATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAA \
        -A GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTTAGTCCACGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAA \
        -o {output.trimmed1} -p {output.trimmed2} \
        {input.reads1} {input.reads2}"
