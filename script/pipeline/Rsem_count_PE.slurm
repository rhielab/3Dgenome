#!/bin/bash
#SBATCH --array=0-2
#SBATCH --job-name=c42b_Rsem.PE.count
#SBATCH --output=c42b_Rsem.PE.count.out
#SBATCH --error=c42b_Rsem.PE.count.err
#SBATCH --time=23:59:59
#SBATCH --mem-per-cpu=16GB
#SBATCH --mail-user=beoungle@usc.edu
#SBATCH --mail-type=ALL

module load perl r bowtie2 star
module load python/3.7.6
module load rsem 

export R_LIBS_USER=/project/rhie_130/suhn/shared/bioinformatic_tools/R4.0_lib
Type='PE'
mkdir /scratch/beoungle/probing_rnaseq/c42b/output/c42b/Rsem_Quantification_PE
cd /scratch/beoungle/probing_rnaseq/c42b/output/c42b/STAR_Alignment_PE
FILES=($(ls *Aligned.toTranscriptome.out.bam))
INPUT=${FILES[$SLURM_ARRAY_TASK_ID]}
Prefix=($(echo ${FILES[$SLURM_ARRAY_TASK_ID]} | rev | cut -d . -f 4- | rev))

if [ $Type == 'SE' ]; then 
rsem-calculate-expression -p 6 --bam --no-bam-output $INPUT /project/rhie_130/suhn/shared/genomes_new/hg38/rsem_index/index /scratch/beoungle/probing_rnaseq/c42b/output/c42b/Rsem_Quantification_PE/$Prefix
else
rsem-calculate-expression -p 6 --bam --no-bam-output $INPUT --paired-end /project/rhie_130/suhn/shared/genomes_new/hg38/rsem_index/index /scratch/beoungle/probing_rnaseq/c42b/output/c42b/Rsem_Quantification_PE/$Prefix
fi