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

computeMatrix reference-point --referencePoint center  -b 1000 -a 1000 -S /project/rhie_130/suhn/beoungle/combined/project/rhie_130/suhn/beoungle/combined/output/mapped_C42B_Novagene_all.PT.bigwig -R CTCF_motif_peak_in_Micro_C_3_Bil.bed -p 16 -o Micro_C_3_Billion.CTCF_motif_in.1kb.RPKM.gz

plotHeatmap -m Micro_C_3_Billion.CTCF_motif_in.1kb.RPKM.gz -o Micro_C_3_Billion.CTCF_motif_in.pdf --legendLocation none --yMin 3 --yMax 5

computeMatrix reference-point --referencePoint center  -b 1000 -a 1000 -S /project/rhie_130/suhn/beoungle/combined/mapped_C42B_Novagene_all.PT.bigwig -R NDR_no_feature_in_Micro_C_3_Bil.chr.bed -p 16 -o NDR_no_feature_in_Micro_C_3_Bil.chr.1kb.RPKM.gz

plotHeatmap -m NDR_no_feature_in_Micro_C_3_Bil.chr.1kb.RPKM.gz -o NDR_no_feature_in_Micro_C_3_Bil.chr.pdf --legendLocation none --yMin 3.8 --yMax 4.8


computeMatrix reference-point --referencePoint center  -b 1000 -a 1000 -S /project/rhie_130/suhn/beoungle/combined/mapped_C42B_Novagene_all.PT.bigwig -R TSS_in_Micro_C_3_Bil.chr.bed -p 16 -o TSS_in_Micro_C_3_Bil.chr.RPKM.gz

plotHeatmap -m TSS_in_Micro_C_3_Bil.chr.RPKM.gz -o TSS_in_Micro_C_3_Bil.chr.pdf --legendLocation none --yMin 3.4 --yMax 4.5



