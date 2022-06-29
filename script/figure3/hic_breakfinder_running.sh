#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=120GB
#SBATCH --time=47:59:59
#SBATCH --mail-user=zexunwu@usc.edu
#SBATCH --mail-type=all
#SBATCH --account=rhie_130
#SBATCH --partition=largemem

hic_breakfinder=/project/rhie_130/suhn/beoungle/neoloop/neoloop_hic/hic_breakfinder/src/hic_breakfinder

#input files
input_bam=/project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/figure6/hic_breakfinder/input/mapped_uni1945.PT.bam
inter_TXT=/project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/figure6/hic_breakfinder/input/inter_expect_1Mb.hg38.txt
intra_TXT=/project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/figure6/hic_breakfinder/input/intra_expect_100kb.hg38.txt
output_path=/project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/figure6/hic_breakfinder/output

cd $output_path
$hic_breakfinder --bam-file $input_bam --exp-file-inter $inter_TXT --exp-file-intra $intra_TXT --name C42B_Capture_Micro_C
#$hic_breakfinder --bam-file mapped_C42B.PT.bam --exp-file-inter inter_expect_1Mb.hg38.txt --exp-file-intra intra_expect_100kb.hg38.txt --name C42B_HiC
#$hic_breakfinder --bam-file mapped_C42B_Novagene_all.PT.bam --exp-file-inter inter_expect_1Mb.hg38.txt --exp-file-intra intra_expect_100kb.hg38.txt --name C42B_3_Billion
#$hic_breakfinder --bam-file mapped_C42B_Micro_2_Billion.PT.bam --exp-file-inter inter_expect_1Mb.hg38.txt --exp-file-intra intra_expect_100kb.hg38.txt --name C42B_2_Billion
