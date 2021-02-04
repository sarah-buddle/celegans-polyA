# local machine (run from workflow folder)
# samtools
conda create --name samtools=1.11 --channel conda-forge --channel bioconda --yes samtools=1.11
conda activate samtools=1.11
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/samtools=1.11.yaml
conda deactivate

# deseq2
conda create --name bioconductor-deseq2=1.30.0 --channel conda-forge --channel bioconda --yes bioconductor-deseq2=1.30.0
conda activate bioconductor-deseq2=1.30.0
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/bioconductor-deseq2=1.30.0.yaml
conda deactivate

# deseq2 and ggplot2
conda create --name bioconductor-deseq2=1.30.0_r-ggplot2=3.3.1 --channel conda-forge --channel bioconda --yes bioconductor-deseq2=1.30.0 r-ggplot2=3.3.1
conda activate bioconductor-deseq2=1.30.0_r-ggplot2=3.3.1
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/bioconductor-deseq2=1.30.0_r-ggplot2=3.3.1.yaml
conda deactivate

# gff3 sort
conda create --name gff3sort=0.1.a1a2bc9 --channel conda-forge --channel bioconda --yes gff3sort=0.1.a1a2bc9
conda activate gff3sort=0.1.a1a2bc9

# rtracklayer
conda create --name bioconductor-rtracklayer=1.50.0 --channel conda-forge --channel bioconda --yes bioconductor-rtracklayer=1.50.0
conda activate bioconductor-rtracklayer=1.50.0
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/bioconductor-rtracklayer=1.50.0.yaml
conda deactivate