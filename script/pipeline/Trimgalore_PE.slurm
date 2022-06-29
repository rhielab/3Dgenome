#!/bin/bash
#SBATCH --array=0-2
#SBATCH --job-name=c42b_TrimGalore.PE
#SBATCH --output=c42b_TrimGalore.PE.out
#SBATCH --error=c42b_TrimGalore.PE.err
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=16GB
#SBATCH --mail-user=beoungle@usc.edu
#SBATCH --mail-type=ALL
module load python/3.7.6
module load fastqc
module load trimgalore 
cd /scratch/beoungle/probing_rnaseq/c42b/output/c42b/c42b_PE_linked_folder_1
FILES_1=($(find `pwd` -type l -name '*gz' | sort)) 
FILES=($(ls *fq *fastq *gz | sort))
cd /scratch/beoungle/probing_rnaseq/c42b/output/c42b/c42b_PE_linked_folder_2
FILES_2=($(find `pwd` -type l -name '*gz' | sort)) 

INPUT_1=${FILES_1[$SLURM_ARRAY_TASK_ID]}
INPUT_2=${FILES_2[$SLURM_ARRAY_TASK_ID]}
OUTPUTDIR=($(echo ${FILES[$SLURM_ARRAY_TASK_ID]} | rev | cut -d . -f 3- | rev))

mkdir /scratch/beoungle/probing_rnaseq/c42b/output/c42b/Trim_PE 

trim_galore --phred33 --fastqc --output_dir /scratch/beoungle/probing_rnaseq/c42b/output/c42b/Trim_PE/$OUTPUTDIR --paired $INPUT_1 $INPUT_2