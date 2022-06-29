#!/bin/bash -e
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=32GB
#SBATCH --time=47:59:59
#SBATCH --mail-user=zexunwu@usc.edu
#SBATCH --mail-type=all
#SBATCH --output=C42B_Micro_1_Billion-4dn_matrices.out
#SBATCH --error=C42B_Micro_1_Billion-4dn_matrices.err
#SBATCH --partition=main

module purge
module load gcc/8.3.0
module load openjdk/11.0.2
module load anaconda3/2021.05; export PYTHONNOUSERSITE=1
eval "$(conda shell.bash hook)"
conda activate /project/rhie_130/suhn/shared/conda_env/4dn_pipeline

cd /project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/4DN_pipeline_testing/output

java -Xmx128000m -Djava.awt.headless=true -jar /project/rhie_130/suhn/shared/juicer/scripts/juicer_tools.jar pre -r 1000  mapped_C42B_Micro_1_Billion.pairs C42B_Micro_1_Billion-1kb.hic /project/rhie_130/suhn/shared/genomes_new/hg38/chrom.sizes
java -Xmx128000m -Djava.awt.headless=true -jar /project/rhie_130/suhn/shared/juicer/scripts/juicer_tools.jar pre  mapped_C42B_Micro_1_Billion.pairs C42B_Micro_1_Billion.hic /project/rhie_130/suhn/shared/genomes_new/hg38/chrom.sizes

bgzip -l 9 --threads 4 -c mapped_C42B_Micro_1_Billion.pairs > mapped_C42B_Micro_1_Billion.pairs.gz
pairix mapped_C42B_Micro_1_Billion.pairs.gz

cooler cload pairix -p 4 /project/rhie_130/suhn/shared/genomes_new/hg38/chrom.sizes:1000 mapped_C42B_Micro_1_Billion.pairs.gz C42B_Micro_1_Billion.cool
cooler zoomify --balance --balance-args "--max-iters 1000 -p 4" -p 4 C42B_Micro_1_Billion.cool

rm -rf "/scratch2/zexunwu/4dn_tmp_701590767785"

