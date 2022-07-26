#To run the script, require the cooltools conda environment. 
module load gcc/8.3.0
eval "$(conda shell.bash hook)"
conda activate /project/rhie_130/suhn/shared/conda_env/cooltools
export PYTHONPATH="${PYTHONPATH}:/home1/zexunwu/.local/lib/python3.8/site-packages"

#1st - mcool file
#2nd - output path
#3rd - prefix for the output
python cooltools_heatmap.py \
input/uni1945_2.5kb_1_billion.filt.mcool \
output/ \
20kb_chr7_31M-43M_Micro_C_1_Bil
