#!/bin/bash
#SBATCH --job-name=c42b_STAR.SE
#SBATCH --output=c42b_STAR.SE.out
#SBATCH --error=c42b_STAR.SE.err
#SBATCH --time=48:00:00
#SBATCH --mem=32GB
#SBATCH --mail-user=beoungle@usc.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
module load star 
mkdir /scratch/beoungle/probing_rnaseq/c42b/output/c42b/STAR_Alignment_SE 
cd /scratch/beoungle/probing_rnaseq/c42b/output/c42b/Trim_SE 
for i in $(ls);
do
cd ${i}
STAR --runThreadN 4 --genomeDir /project/rhie_130/suhn/shared/genomes_new/hg38/star_index --readFilesIn *.fq.gz --outFileNamePrefix /scratch/beoungle/probing_rnaseq/c42b/output/c42b/STAR_Alignment_SE/${i} --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --quantMode TranscriptomeSAM GeneCounts --readFilesCommand zcat
cd ..
done