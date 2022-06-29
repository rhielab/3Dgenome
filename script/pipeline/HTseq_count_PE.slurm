#!/bin/bash
#SBATCH --array=0-2
#SBATCH --job-name=c42b_HTseq.PE.count
#SBATCH --output=c42b_HTseq.count.PE.out
#SBATCH --error=c42b_HTseq.count.PE.err
#SBATCH --time=23:59:59
#SBATCH --mem-per-cpu=16GB
#SBATCH --mail-user=beoungle@usc.edu
#SBATCH --mail-type=ALL

module load python/3.7.6
cd /scratch/beoungle/probing_rnaseq/c42b/output/c42b
mkdir HTseq_quantification
cd /scratch/beoungle/probing_rnaseq/c42b/output/c42b/STAR_Alignment_PE
FILES=($(ls *.sortedByCoord.out.bam))
INPUT=${FILES[$SLURM_ARRAY_TASK_ID]}
Prefix=($(echo ${FILES[$SLURM_ARRAY_TASK_ID]} | rev | cut -d . -f 4- | rev))

htseq-count -s no -r pos -f bam $INPUT /project/rhie_130/suhn/shared/genomes_new/hg38/annotations/gencode.gtf > /scratch/beoungle/probing_rnaseq/c42b/output/c42b/HTseq_quantification/${Prefix}.counts