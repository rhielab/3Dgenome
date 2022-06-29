#!/bin/bash
#SBATCH --job-name=c42b_GenesCount
#SBATCH --output=c42b_GenesCount.out
#SBATCH --error=c42b_GenesCount.err
#SBATCH --time=23:59:59
#SBATCH --mem-per-cpu=16GB
#SBATCH --mail-user=beoungle@usc.edu
#SBATCH --mail-type=ALL

module load r
export R_LIBS_USER=/project/rhie_130/suhn/shared/bioinformatic_tools/R4.0_lib
Rscript /scratch/beoungle/probing_rnaseq/c42b/c42b_RNA-seq_scripts/GeneCount.R