#!/bin/bash -e
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=32GB
#SBATCH --time=23:59:59
#SBATCH --mail-user=zexunwu@discovery.usc.edu
#SBATCH --mail-type=all
#SBATCH --output=uni1945_2_billion-TopDom.out
#SBATCH --error=uni1945_2_billion-TopDom.err
#SBATCH --partition=main

module purge
module load gcc/8.3.0
module load anaconda3/2021.05; export PYTHONNOUSERSITE=1
eval "$(conda shell.bash hook)"
conda activate /project/rhie_130/suhn/shared/conda_env/tad_interaction_pipeline

module load usc
module load r/4.0.3
export R_LIBS_USER=/project/rhie_130/suhn/shared/bioinformatic_tools/R4.0_lib

#Generate sparse matrix format using python script
python uni1945_2_billion-straw.py /project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/TAD_pipeline_test/input/C42B_Novagene_all.hic 50000 uni1945_2_billion

#Generate topdom output
Rscript uni1945_2_billion-topdom.R uni1945_2_billion 50000

#Combine topdom output and produce final output topdom bed file with domain information
tail -n +2 *topdom.bed | grep -v '==>' > uni1945_2_billion-50000-topdom.bed



