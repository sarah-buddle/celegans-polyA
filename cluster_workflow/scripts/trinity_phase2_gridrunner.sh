#!/bin/bash -e
#SBATCH --time=30:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --hint=nomultithread
#SBATCH -p IACT
#SBATCH --job-name=trinity_phase2_gridrunner_bristol_pf_rep1
#SBATCH -o polyA/reference_free/trinity/bristol/pf/rep1/trinity_output_phase2
#SBATCH -e polyA/reference_free/trinity/bristol/pf/rep1/trinity_error_phase2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sb2226@cam.ac.uk
srun Trinity --CPU ${SLURM_CPUS_PER_TASK} --max_memory 20G \
--grid_exec "${SLURM_SUBMIT_DIR}/HpcGridRunner-1.0.2/hpc_cmds_GridRunner.pl --grid_conf ${SLURM_SUBMIT_DIR}/HpcGridRunner-1.0.2/grid_config.conf -c" \
--seqType fq \
--left polyA/QC/trimmed_fastq/bristol/pf/rep1/reads1/bristol_pf_rep1_1_trimmed.fastq \
--right polyA/QC/trimmed_fastq/bristol/pf/rep1/reads2/bristol_pf_rep1_2_trimmed.fastq \
--output polyA/reference_free/trinity/bristol/pf/rep1/trinity
