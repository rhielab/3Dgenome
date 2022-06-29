#!/bin/bash -e
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64GB
#SBATCH --time=23:59:59
#SBATCH --mail-user=zexunwu@discovery.usc.edu
#SBATCH --mail-type=all
#SBATCH --output=uni1945_2_billion-Mustache.out
#SBATCH --error=uni1945_2_billion-Mustache.err
#SBATCH --partition=main


module purge
module load gcc/8.3.0
module load anaconda3/2021.05; export PYTHONNOUSERSITE=1
eval "$(conda shell.bash hook)"
conda activate /project/rhie_130/suhn/shared/conda_env/tad_interaction_pipeline 

#Check the name of chromosome/genome whether it starts with chr or just the number
chr_prefix=

#Generate mustache loops for each chromosomes, then combine at the end. Check whether resolution is less than 5kb or less, and if less than 5kb, then use specialized setting for higher resolution loops.
for resolution in $(echo 50kb,10kb | tr , " ")
do
chr_test=$(/project/rhie_130/suhn/shared/bioinformatic_tools/straw NONE /project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/TAD_pipeline_test/input/C42B_Novagene_all.hic 1:1:1 1:1:1 BP $resolution 2>&1 /dev/null || true)
if [ "$chr_test" == "One of the chromosomes wasn't found in the file. Check that the chromosome name matches the genome." ];
then
chr_prefix="chr"
fi
for chr_number in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
do
test_resolution=$(echo $resolution | tr -d kb)
if (( $test_resolution < 5 ));
then
mustache -f /project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/TAD_pipeline_test/input/C42B_Novagene_all.hic -ch $chr_prefix$chr_number -r $resolution -st 0.7 -pt 0.1 -o uni1945_2_billion_Mustache-$resolution-chr$chr_number-loop.tsv
else
mustache -f /project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/TAD_pipeline_test/input/C42B_Novagene_all.hic -ch $chr_prefix$chr_number -r $resolution -o uni1945_2_billion_Mustache-$resolution-chr$chr_number-loop.tsv
fi
done
mustache -f /project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/TAD_pipeline_test/input/C42B_Novagene_all.hic -ch ${chr_prefix}Y -r $resolution -o uni1945_2_billion_Mustache-$resolution-chrY-loop.tsv || echo "Chromosome Y Failed"
tail -n +2 uni1945_2_billion_Mustache-$resolution-chr*-loop.tsv | grep '000' > uni1945_2_billion_Mustache-$resolution-all-loop.tsv
if [ $(wc -l < uni1945_2_billion_Mustache-$resolution-all-loop.tsv) -le 1 ];
then
echo "Not enough reads in the file for this resolution"
exit 1
fi
done


