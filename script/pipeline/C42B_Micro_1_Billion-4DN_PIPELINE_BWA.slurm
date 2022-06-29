#!/bin/bash -e
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=64GB
#SBATCH --time=167:59:59
#SBATCH --mail-user=zexunwu@usc.edu
#SBATCH --mail-type=all
#SBATCH --output=C42B_Micro_1_Billion-4dn_bwa.out
#SBATCH --error=C42B_Micro_1_Billion-4dn_bwa.err
#SBATCH --partition=rhie_130

mkdir -p /scratch2/zexunwu/4dn_tmp_701590767785/samtools;
mkdir -p /scratch2/zexunwu/4dn_tmp_701590767785/pairtools;

module purge
module load gcc/8.3.0
module load anaconda3/2021.05; export PYTHONNOUSERSITE=1
eval "$(conda shell.bash hook)"
conda activate /project/rhie_130/suhn/shared/conda_env/4dn_pipeline 

cd /project/rhie_130/suhn/beoungle/pipeline_dovetail/c42b_1_billion
shopt -s nullglob # avoid error for nonexistent patterns
FILES_1=($(ls {*{R,_}1.fastq.gz,*{R,_}1.fastq,*{R,_}1.fq,*{R,_}1.fq.gz}))
FILES_2=($(ls {*{R,_}2.fastq.gz,*{R,_}2.fastq,*{R,_}2.fq,*{R,_}2.fq.gz}))
shopt -u nullglob

if [ ${#FILES_1[@]} -eq 0 ]; then
	echo "ERROR: No R1 fastqs found!"
	exit 1
fi
if [ ${#FILES_2[@]} -eq 0 ]; then
	echo "ERROR: No R2 fastqs found!"
	exit 1
fi

mkdir -p /project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/4DN_pipeline_testing/output

bwa mem -5SP -T0 -t4 /project/rhie_130/suhn/shared/genomes_new/hg38/bwa_index/index <(cat ${FILES_1[@]:0}) <(cat ${FILES_2[@]:0})  | samtools view -bS - > /project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/4DN_pipeline_testing/output/aligned_C42B_Micro_1_Billion.bam
