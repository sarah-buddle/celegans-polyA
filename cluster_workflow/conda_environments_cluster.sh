# Create all necessary conda environments for cluster workflow

# samtools
conda create --name samtools=1.11 --channel conda-forge --channel bioconda --yes samtools=1.11
conda activate samtools=1.11
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/samtools=1.11.yaml
conda deactivate

# fastq-tools
conda create --name fastq-tools=0.8.3 --channel conda-forge --channel bioconda --yes fastq-tools=0.8.3
conda activate fastq-tools=0.8.3
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/fastq-tools=0.8.3.yaml
conda deactivate

# seqtk
conda create --name seqtk=1.3 --channel conda-forge --channel bioconda --yes seqtk=1.3
conda activate seqtk=1.3
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/seqtk=1.3.yaml
conda deactivate

# fastqc
conda create --name fastqc=0.11.9 --channel conda-forge --channel bioconda --yes fastqc=0.11.9
conda activate fastqc=0.11.9
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/fastqc=0.11.9.yaml
conda deactivate

# multiqc
conda create --name multiqc=1.9 --channel conda-forge --channel bioconda --yes multiqc=1.9
conda activate multiqc=1.9
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/multiqc=1.9.yaml
conda deactivate

# cutadapt
conda create --name cutadapt=3.0 --channel conda-forge --channel bioconda --yes cutadapt=3.0
conda activate cutadapt=3.0
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/cutadapt=3.0.yaml
conda deactivate

# Trinity
conda create --name Trinity=2.11.0 --channel conda-forge --channel bioconda --yes Trinity=2.11.0
conda activate Trinity=2.11.0
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/Trinity=2.11.0.yaml
conda deactivate

# hisat2
conda create --name hisat2=2.2.1 --channel conda-forge --channel bioconda --yes hisat2=2.2.1
conda activate hisat2=2.2.1
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/hisat2=2.2.1.yaml
conda deactivate

# stringtie
conda create --name stringtie=2.1.4 --channel conda-forge --channel bioconda --yes stringtie=2.1.4
conda activate stringtie=2.1.4
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/stringtie=2.1.4.yaml
conda deactivate

# maker
conda create --name maker=2.31.10 --channel conda-forge --channel bioconda --yes maker=2.31.10
conda activate maker=2.31.10
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/maker=2.31.10.yaml
conda deactivate

# repeatmasker
conda create --name repeatmasker=4.1.1 --channel conda-forge --channel bioconda --yes repeatmasker=4.1.1
conda activate repeatmasker=4.1.1
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/repeatmasker=4.1.1.yaml
conda deactivate

# htseq-count
conda create --name htseq=0.12.4 --channel conda-forge --channel bioconda --yes htseq=0.12.4
conda activate htseq=0.12.4
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/htseq=0.12.4.yaml
conda deactivate

# gffcompare
conda create --name gffcompare=0.11.2 --channel conda-forge --channel bioconda --yes gffcompare=0.11.2
conda activate gffcompare=0.11.2
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/gffcompare=0.11.2.yaml
conda deactivate

# soapdenovo-trans
conda create --name soapdenovo-trans=1.04 --channel conda-forge --channel bioconda --yes soapdenovo-trans=1.04
conda activate soapdenovo-trans=1.04
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/soapdenovo-trans=1.04.yaml
conda deactivate
