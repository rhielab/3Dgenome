#!/bin/bash -e
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=32GB
#SBATCH --time=47:59:59
#SBATCH --mail-user=zexunwu@usc.edu
#SBATCH --mail-type=all
#SBATCH --output=C42B_Micro_1_Billion-samtools.out
#SBATCH --error=C42B_Micro_1_Billion-samtools.err
#SBATCH --partition=main

module purge
module load gcc/8.3.0
module load anaconda3/2021.05; export PYTHONNOUSERSITE=1
eval "$(conda shell.bash hook)"
conda activate /project/rhie_130/suhn/shared/conda_env/4dn_pipeline

cd /project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/4DN_pipeline_testing/output

#Generate final bam file
samtools sort -@ 4 -m 16G -T /scratch2/zexunwu/4dn_tmp_701590767785/samtools -o mapped_C42B_Micro_1_Billion.PT.bam unsorted_C42B_Micro_1_Billion.bam

#Index for future use
samtools index -@ 3 mapped_C42B_Micro_1_Billion.PT.bam # in index, for some reason, -@ is the number of ADDITIONAL threads

conda deactivate

conda activate /project/rhie_130/suhn/shared/conda_env/deeptools

#Create bigwig file
bamCoverage -p 4 --normalizeUsing RPKM --ignoreDuplicates -e 0 -bs 10 -b mapped_C42B_Micro_1_Billion.PT.bam -o mapped_C42B_Micro_1_Billion.PT.bigwig

rm unsorted_C42B_Micro_1_Billion.bam
