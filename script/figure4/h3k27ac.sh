#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8GB
#SBATCH --time=23:59:59
#SBATCH --mail-user=beoungle@usc.edu
#SBATCH --mail-type=all#!/bin/bash

module load gcc/8.3.0
module load anaconda3

eval "$(conda shell.bash hook)"

conda activate deeptools

TMPDIR=/project/rhie_130/suhn/beoungle/combined/temp

computeMatrix reference-point --referencePoint center  -b 1000 -a 1000 -S /project/rhie_130/suhn/beoungle/combined/mapped_C42B_Novagene_all.PT.bigwig -R H3K27ac_peak_in_Micro_C_3_Bil.chrndr.bed -p 16 -o H3K27ac_peak_in_Micro_C_3_Bil.chrndr.1kb.RPKM.gz

plotHeatmap -m H3K27ac_peak_in_Micro_C_3_Bil.chrndr.1kb.RPKM.gz -o H3K27ac_peak_in_Micro_C_3_Bil.chrndr.pdf --legendLocation none --yMin 3.25 --yMax 4.75

