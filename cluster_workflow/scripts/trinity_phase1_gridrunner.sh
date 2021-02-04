#!/bin/bash -e
#SBATCH --time=30:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --hint=nomultithread
#SBATCH -p IACT
#SBATCH --job-name=trinity_phase1_bristol_pf_rep1
#SBATCH -o polyA/reference_free/trinity/bristol/pf/rep1/trinity_output_phase1
#SBATCH -e polyA/reference_free/trinity/bristol/pf/rep1/trinity_error_phase1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sb2226@cam.ac.uk

srun Trinity --no_distributed_trinity_exec --CPU ${SLURM_CPUS_PER_TASK} --max_memory 100G \
--seqType fq \
--left polyA/QC/trimmed_fastq/bristol/pf/rep1/reads1/bristol_pf_rep1_1_trimmed.fastq \
--right polyA/QC/trimmed_fastq/bristol/pf/rep1/reads2/bristol_pf_rep1_2_trimmed.fastq \
--output polyA/reference_free/trinity/bristol/pf/rep1/trinity
