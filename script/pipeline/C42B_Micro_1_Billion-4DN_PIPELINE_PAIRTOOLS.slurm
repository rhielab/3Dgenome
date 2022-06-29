#!/bin/bash -e
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=32GB
#SBATCH --time=47:59:59
#SBATCH --mail-user=zexunwu@usc.edu
#SBATCH --mail-type=all
#SBATCH --output=C42B_Micro_1_Billion-4dn_pairtools.out
#SBATCH --error=C42B_Micro_1_Billion-4dn_pairtools.err
#SBATCH --partition=main

module purge
module load gcc/8.3.0
module load anaconda3/2021.05; export PYTHONNOUSERSITE=1
eval "$(conda shell.bash hook)"
conda activate /project/rhie_130/suhn/shared/conda_env/4dn_pipeline

cd /project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/4DN_pipeline_testing/output

#record valid ligation events
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 4 --nproc-out 4 --chroms-path /project/rhie_130/suhn/shared/genomes_new/hg38/chrom.sizes aligned_C42B_Micro_1_Billion.bam > parsed_C42B_Micro_1_Billion.pairsam

#sort the pairsam file
pairtools sort --nproc 4 --tmpdir=/scratch2/zexunwu/4dn_tmp_701590767785/pairtools parsed_C42B_Micro_1_Billion.pairsam > sorted_C42B_Micro_1_Billion.pairsam

#remove PCR duplicates
pairtools dedup --nproc-in 4 --nproc-out 4 --mark-dups --output-stats C42B_Micro_1_Billion-stats.txt --output dedup_C42B_Micro_1_Billion.pairsam sorted_C42B_Micro_1_Billion.pairsam

#Generate summarized stats file (Added 20220421 Ethan)
python3 /project/rhie_130/suhn/shared/4DN.Pipeline/get_qc.py -p C42B_Micro_1_Billion-stats.txt > C42B_Micro_1_Billion-summarizedStats.txt

#Generate .pair and bam files
pairtools split --nproc-in 4 --nproc-out 4 --output-pairs mapped_C42B_Micro_1_Billion.pairs --output-sam unsorted_C42B_Micro_1_Billion.bam dedup_C42B_Micro_1_Billion.pairsam

#Generate HiC-Pro validpair file
grep -v '#' mapped_C42B_Micro_1_Billion.pairs | awk -F "	" '{print $1"	"$2"	"$3"	"$6"	"$4"	"$5"	"$7}' > C42B_Micro_1_Billion.validPairs

rm parsed_C42B_Micro_1_Billion.pairsam sorted_C42B_Micro_1_Billion.pairsam dedup_C42B_Micro_1_Billion.pairsam
