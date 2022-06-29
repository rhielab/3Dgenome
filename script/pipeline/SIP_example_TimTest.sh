#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=64GB
#SBATCH --time=71:59:59
#SBATCH --mail-user=zexunwu@usc.edu
#SBATCH --mail-type=all
#SBATCH --partition=rhie
#SBATCH --account=rhie_130

module load openjdk


#SIP also uses java, and juicertools as well.
#Each resolution needs to be ran separately.

#java -jar /project/rhie_130/suhn/shared/bioinformatic_tools/SIP/SIP_HiC_v1.6.1.jar hic /project/rhie_130/suhn/beoungle/combined/project/rhie_130/suhn/beoungle/combined/output/C42B_Novagene_all-1kb.hic /project/rhie_130/suhn/shared/genomes_new/hg38/hg38.chrom.sizes /project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/SIP_test/output/C42B_Micro-C_3_Billion_1kb /project/rhie_130/suhn/shared/juicer/scripts/juicer_tools.jar -res 1000 -fdr 0.10

java -jar /project/rhie_130/suhn/shared/bioinformatic_tools/SIP/SIP_HiC_v1.6.1.jar hic /project/rhie_130/suhn/beoungle/combined/project/rhie_130/suhn/beoungle/combined/output/C42B_Novagene_all.hic /project/rhie_130/suhn/shared/genomes_new/hg38/hg38.chrom.sizes /project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/SIP_test/output/C42B_Micro-C_3_Billion_50kb /project/rhie_130/suhn/shared/juicer/scripts/juicer_tools.jar -res 50000 -fdr 0.10

java -jar /project/rhie_130/suhn/shared/bioinformatic_tools/SIP/SIP_HiC_v1.6.1.jar hic /project/rhie_130/suhn/beoungle/combined/project/rhie_130/suhn/beoungle/combined/output/C42B_Novagene_all.hic /project/rhie_130/suhn/shared/genomes_new/hg38/hg38.chrom.sizes /project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/SIP_test/output/C42B_Micro-C_3_Billion_25kb /project/rhie_130/suhn/shared/juicer/scripts/juicer_tools.jar -res 25000 -fdr 0.10

java -jar /project/rhie_130/suhn/shared/bioinformatic_tools/SIP/SIP_HiC_v1.6.1.jar hic /project/rhie_130/suhn/beoungle/combined/project/rhie_130/suhn/beoungle/combined/output/C42B_Novagene_all.hic /project/rhie_130/suhn/shared/genomes_new/hg38/hg38.chrom.sizes /project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/SIP_test/output/C42B_Micro-C_3_Billion_10kb /project/rhie_130/suhn/shared/juicer/scripts/juicer_tools.jar -res 10000 -fdr 0.10

java -jar /project/rhie_130/suhn/shared/bioinformatic_tools/SIP/SIP_HiC_v1.6.1.jar hic /project/rhie_130/suhn/beoungle/combined/project/rhie_130/suhn/beoungle/combined/output/C42B_Novagene_all.hic /project/rhie_130/suhn/shared/genomes_new/hg38/hg38.chrom.sizes /project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/SIP_test/output/C42B_Micro-C_3_Billion_5kb /project/rhie_130/suhn/shared/juicer/scripts/juicer_tools.jar -res 5000 -fdr 0.10

java -jar /project/rhie_130/suhn/shared/bioinformatic_tools/SIP/SIP_HiC_v1.6.1.jar hic /project/rhie_130/suhn/beoungle/combined/project/rhie_130/suhn/beoungle/combined/output/C42B_Novagene_all-2kb.hic /project/rhie_130/suhn/shared/genomes_new/hg38/hg38.chrom.sizes /project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/SIP_test/output/C42B_Micro-C_3_Billion_2kb /project/rhie_130/suhn/shared/juicer/scripts/juicer_tools.jar -res 2000 -fdr 0.10

java -jar /project/rhie_130/suhn/shared/bioinformatic_tools/SIP/SIP_HiC_v1.6.1.jar hic /project/rhie_130/suhn/beoungle/combined/project/rhie_130/suhn/beoungle/combined/output/C42B_Novagene_all-1kb.hic /project/rhie_130/suhn/shared/genomes_new/hg38/hg38.chrom.sizes /project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/SIP_test/output/C42B_Micro-C_3_Billion_1kb /project/rhie_130/suhn/shared/juicer/scripts/juicer_tools.jar -res 1000 -fdr 0.10
