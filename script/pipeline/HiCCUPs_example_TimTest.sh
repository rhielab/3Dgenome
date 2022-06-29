#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32GB
#SBATCH --time=47:59:59
#SBATCH --mail-user=zexunwu@usc.edu
#SBATCH --mail-type=all
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100:1

module load openjdk
module load jdk

#GPU-enabled software often requires the CUDA Toolkit or the cuDNN library.
module load gcc/8.3.0
module load cuda/10.1.243
module load cudnn/8.0.2-10.1

#This is for HiCCUPs, chromatin interaction program built with Juicer.
cd /project/rhie_130/suhn/shared/juicer/scripts

#HiCCUPs uses juicertools; hiccups is a tool contained with juicertools.
#I had to lower some setting for 2kb and 1kb, because with normal setting, I couldn't identify interactions at all with hiccups. Notably, you can look at the parameters that I specified between 5000,10000,25000 to notice the differences in the setting. -f is the only exception.
#HiCCUPS didn't work with 50kb so we don't have 50,000 for -r.
#java -jar juicer_tools.jar hiccups -r 5000,10000,25000 -f 0.2,0.2,0.2 /project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/HiCCUP_test/input/C42B_Novagene_all.hic /project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/HiCCUP_test/output
java -jar juicer_tools.jar hiccups --threads 1 -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y -r 1000 -f 0.2 -p 20 -i 50 -d 20000 /project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/HiCCUP_test/input/C42B_Novagene_all-1kb.hic /project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/HiCCUP_test/output
java -jar juicer_tools.jar hiccups --threads 1 -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y -r 2000 -f 0.2 -p 10 -i 25 -d 20000 /project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/HiCCUP_test/input/C42B_Novagene_all-2kb.hic /project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/HiCCUP_test/output
