#!/bin/bash
#SBATCH -n 24
#SBATCH -N 1
#SBATCH -p IACT
#SBATCH --job-name=trinity_bristol_as_rep1
#SBATCH -o polyA/reference_free/trinity/bristol/as/rep1/trinity_output
#SBATCH -e polyA/reference_free/trinity/bristol/as/rep1/trinity_error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sb2226@cam.ac.uk
Trinity --seqType fq --max_memory 60G \
--left polyA/QC/trimmed_fastq/bristol/as/rep1/reads1/bristol_as_rep1_1_trimmed.fastq \
--right polyA/QC/trimmed_fastq/bristol/as/rep1/reads2/bristol_as_rep1_2_trimmed.fastq \
--output polyA/reference_free/trinity/bristol/as/rep1/trinity \
--CPU 24
