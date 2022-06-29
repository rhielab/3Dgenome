#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16GB
#SBATCH --time=47:59:59
#SBATCH --mail-user=zexunwu@usc.edu
#SBATCH --mail-type=all
#SBATCH --partition=main
#SBATCH --job-name=NeoLoop_test

module load gcc/8.3.0
module load anaconda3
eval "$(conda shell.bash hook)"
conda activate /project/rhie_130/suhn/shared/conda_env/neoloop
#The below file is from BE, I test it to ensure that whether the commands below are changing the original cool file or not
INPUT_COOL=/project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/figure6/neoloopfinder/C42B_Micro_C_3_Bil_5kb_original.cool
INPUT_SV=/project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/figure6/neoloopfinder/input/C42B_3_Billion.sv.txt
OUTPUT_PATH=/project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/figure6/neoloopfinder/output_2

cd $OUTPUT_PATH
cooler balance -p 1 $INPUT_COOL
#calculate-cnv -H $INPUT_COOL -g hg38 -e uniform --output C42B_Micro_3_Billion_5kb --cachefolder test
#segment-cnv --cnv-file C42B_Micro_3_Billion_5kb.bedgraph --binsize 10000 --output C42B_Micro_3_Billion_5kb_segment
#correct-cnv -H $INPUT_COOL --cnv-file C42B_Micro_3_Billion_5kb_segment.bedgraph
assemble-complexSVs -O C42B_Micro_3_Billion_5kb -B $INPUT_SV -H $INPUT_COOL --balance-type ICE
neoloop-caller -O C42B_Micro_3_Billion_5kb.neo-loops.txt -H $INPUT_COOL --assembly C42B_Micro_3_Billion_5kb.assemblies.txt --no-clustering --prob 0.95 --balance-type ICE
 
