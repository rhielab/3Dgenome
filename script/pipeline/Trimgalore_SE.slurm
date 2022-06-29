#!/bin/bash
#SBATCH --array=0-7
#SBATCH --job-name=c42b_TrimGalore.SE
#SBATCH --output=c42b_TrimGalore.SE.out
#SBATCH --error=c42b_TrimGalore.SE.err
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=16GB
#SBATCH --mail-user=beoungle@usc.edu
#SBATCH --mail-type=ALL
module load python/3.7.6
module load fastqc
module load trimgalore 
cd /scratch/beoungle/probing_rnaseq/c42b/output/c42b/c42b_SE_linked_folder
FILES=($(ls {*fastq.gz,*.fq.gz})) 
INPUT=${FILES[$SLURM_ARRAY_TASK_ID]} 
OUTPUTDIR=($(echo ${FILES[$SLURM_ARRAY_TASK_ID]} | rev | cut -d . -f 3- | rev)) 

mkdir /scratch/beoungle/probing_rnaseq/c42b/output/c42b/Trim_SE 

trim_galore --phred33 --fastqc --output_dir /scratch/beoungle/probing_rnaseq/c42b/output/c42b/Trim_SE/$OUTPUTDIR $INPUT