# conda environments needed for workflow on local machine

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
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/gff3sort=0.1.a1a2bc9.yaml
conda deactivate

# rtracklayer
conda create --name bioconductor-rtracklayer=1.50.0 --channel conda-forge --channel bioconda --yes bioconductor-rtracklayer=1.50.0
conda activate bioconductor-rtracklayer=1.50.0
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/bioconductor-rtracklayer=1.50.0.yaml
conda deactivate

# tidyverse
conda create --name r-tidyverse=1.2.1 --yes r-tidyverse=1.2.1
conda activate r-tidyverse=1.2.1
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/r-tidyverse=1.2.1.yaml
conda deactivate

# dplyr
conda create --name r-dplyr=0.8.0.1 --yes r-dplyr=0.8.0.1
conda activate r-dplyr=0.8.0.1
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/r-dplyr=0.8.0.1.yaml
conda deactivate

# topgo
conda create --name bioconductor-topgo=2.42.0 --channel conda-forge --channel bioconda --yes bioconductor-topgo=2.42.0
conda activate bioconductor-topgo=2.42.0
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/bioconductor-topgo=2.42.0.yaml
conda deactivate

# GenomicRanges
conda create --name bioconductor-genomicranges=1.42.0 --channel conda-forge --channel bioconda --yes bioconductor-genomicranges=1.42.0
conda activate bioconductor-genomicranges=1.42.0
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/bioconductor-genomicranges=1.42.0.yaml
conda deactivate

# ggplot2
conda create --name r-ggplot2=3.3.1 --channel conda-forge --channel bioconda --yes r-ggplot2=3.3.1
conda activate r-ggplot2=3.3.1
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/r-ggplot2=3.3.1.yaml
conda deactivate

# fqtools
conda create --name fqtools=2.0 --channel conda-forge --channel bioconda --yes fqtools=2.0
conda activate fqtools=2.0

# Genomic Ranges and plyranges
conda create --name bioconductor-genomicranges=1.42.0_bioconductor-plyranges=1.10.0 --channel conda-forge --channel bioconda --yes bioconductor-genomicranges=1.42.0 bioconductor-plyranges=1.10.0
conda activate bioconductor-genomicranges=1.42.0_bioconductor-plyranges=1.10.0
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/bioconductor-genomicranges=1.42.0_bioconductor-plyranges=1.10.0.yaml
conda deactivate

# plyranges
conda create --name bioconductor-plyranges=1.10.0 --channel conda-forge --channel bioconda --yes bioconductor-plyranges=1.10.0
conda activate bioconductor-plyranges=1.10.0
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/bioconductor-plyranges=1.10.0.yaml
conda deactivate

# bedops
conda create --name bedops=2.4.39 --channel conda-forge --channel bioconda --yes bedops=2.4.39
conda activate bedops=2.4.39
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/bedops=2.4.39.yaml
conda deactivate

# igv tools
conda create --name igvtools=2.5.3 --channel conda-forge --channel bioconda --yes igvtools=2.5.3
conda activate igvtools=2.5.3
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/igvtools=2.5.3.yaml
conda deactivate

# blast
conda create --name blast=2.10.1 --channel conda-forge --channel bioconda --yes blast=2.10.1
conda activate blast=2.10.1
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/blast=2.10.1.yaml
conda deactivate

# r-seqinr
conda create --name r-seqinr=3.4_5 --channel conda-forge --channel bioconda --yes r-seqinr=3.4_5
conda activate r-seqinr=3.4_5
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/r-seqinr=3.4_5.yaml
conda deactivate

# rrvgo
conda create --name bioconductor-rrvgo=1.2.0_bioconductor-org.ce.eg.db=3.12.0 \
--channel conda-forge --channel bioconda --yes bioconductor-rrvgo=1.2.0 \
bioconductor-org.ce.eg.db=3.12.0
conda activate bioconductor-rrvgo=1.2.0_bioconductor-org.ce.eg.db=3.12.0
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/bioconductor-rrvgo=1.2.0_bioconductor-org.ce.eg.db=3.12.0.yaml
conda deactivate

# rtracklayer, granges and plyranges
conda create --name bioconductor-rtracklayer=1.50.0_bioconductor-genomicranges=1.42.0_bioconductor-plyranges=1.10.0 \
--channel conda-forge --channel bioconda --yes bioconductor-rtracklayer=1.50.0 \
bioconductor-genomicranges=1.42.0 bioconductor-plyranges=1.10.0
conda activate bioconductor-rtracklayer=1.50.0_bioconductor-genomicranges=1.42.0_bioconductor-plyranges=1.10.0
conda env export --no-builds | sed '$d' | sed '$d' > \
envs/conda/bioconductor-rtracklayer=1.50.0_bioconductor-genomicranges=1.42.0_bioconductor-plyranges=1.10.0.yaml
conda deactivate

# rtracklayer, ggplot and plyranges
conda create --name bioconductor-rtracklayer=1.50.0_r-ggplot2=3.3.1_bioconductor-plyranges=1.10.0 \
--channel conda-forge --channel bioconda --yes bioconductor-rtracklayer=1.50.0 \
r-ggplot2=3.3.1 bioconductor-plyranges=1.10.0
conda activate bioconductor-rtracklayer=1.50.0_r-ggplot2=3.3.1_bioconductor-plyranges=1.10.0
conda env export --no-builds | sed '$d' | sed '$d' > \
envs/conda/bioconductor-rtracklayer=1.50.0_r-ggplot2=3.3.1_bioconductor-plyranges=1.10.0.yaml
conda deactivate

# deseq2 and plyranges
conda create --name bioconductor-deseq2=1.30.0_bioconductor-plyranges=1.10.0 \
--channel conda-forge --channel bioconda --yes bioconductor-deseq2=1.30.0 \
bioconductor-plyranges=1.10.0
conda activate bioconductor-deseq2=1.30.0_bioconductor-plyranges=1.10.0
conda env export --no-builds | sed '$d' | sed '$d' > \
envs/conda/bioconductor-deseq2=1.30.0_bioconductor-plyranges=1.10.0.yaml
conda deactivate

# ggplot2, rggpubr, gridextra
conda create --name r-ggplot2=3.3.1_r-ggpubr=0.4.0_r-gridextra=2.3 \
--channel conda-forge --channel bioconda --yes r-ggplot2=3.3.1 r-ggpubr=0.4.0 r-gridextra=2.3
conda activate r-ggplot2=3.3.1_r-ggpubr=0.4.0_r-gridextra=2.3
conda env export --no-builds | sed '$d' | sed '$d' > \
r-ggplot2=3.3.1_rggpubr=0.4.0_r-gridextra=2.3.yaml
conda deactivate

conda create --name bioconductor-ggtree=2.4.1 \
--channel conda-forge --channel bioconda --yes bioconductor-ggtree=2.4.1
conda activate bioconductor-ggtree=2.4.1
conda env export --no-builds | sed '$d' | sed '$d' > \
envs/conda/bioconductor-ggtree=2.4.1.yaml
conda deactivate
