#!/bin/bash
#SBATCH --job-name=c42b_STAR.PE
#SBATCH --output=c42b_STAR.PE.out
#SBATCH --error=c42b_STAR.PE.err
#SBATCH --time=48:00:00
#SBATCH --mem=32GB
#SBATCH --mail-user=beoungle@usc.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
module load star 
mkdir /scratch/beoungle/probing_rnaseq/c42b/output/c42b/STAR_Alignment_PE 
cd /scratch/beoungle/probing_rnaseq/c42b/output/c42b/Trim_PE 

OUTPUTDIR=($(ls ./ | sort)) 
FILES_1=($(find ./ -type f -name '*1.fq.gz' | sort))
FILES_2=($(find ./ -type f -name '*2.fq.gz' | sort))

for i in {0..2};
do
STAR --runThreadN 4 --genomeDir /project/rhie_130/suhn/shared/genomes_new/hg38/star_index --readFilesIn ${FILES_1[i]} ${FILES_2[i]} --outFileNamePrefix /scratch/beoungle/probing_rnaseq/c42b/output/c42b/STAR_Alignment_PE/${OUTPUTDIR[i]} --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --quantMode TranscriptomeSAM GeneCounts --readFilesCommand zcat
done