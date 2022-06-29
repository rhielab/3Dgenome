#!/bin/bash -e
# $1 out prefix 
# $2 hcg bw 
# $3 gch bw
# $4 ndr bed 
# $5 short name like 'CTCF'
module purge
module load gcc/8.3.0 openblas/0.3.8 r/4.0.3 picard/2.20.8 samtools/1.10
export R_LIBS_USER=/project/rhie_130/suhn/shared/bioinformatic_tools/R4.0_lib
export BISTOOLS=/project/rhie_130/suhn/shared/bioinformatic_tools/Bis-tools_0.90
#export PATH=/project/rhie_130/suhn/shared/bioinformatic_tools/ucsc:$PATH
Dir="`dirname "$(realpath "$1")"`/" # / needed to work around issue
Name="$5" # I think this should be just like RhieBXXX
#Prefix1="`realpath "$1"`" # full path doesnt work
Prefix1="$1"
pushd $Dir
#ln -sf "$4" . # needs to be absolute path in argument
cp "$4" .
popd
L1="`basename $4`"

N1="CTCF"

#awk 'BEGIN { OFS="\t"} {print $1, $2, $3, $4*100}' $2 > $2.percent.bedGraph
#awk 'BEGIN { OFS="\t"} {print $1, $2, $3, $4*100}' $3 > $3.percent.bedGraph

#bedGraphToBigWig $2.percent.bedGraph /project/rhie_130/suhn/shared/genomes_new/$5/sequence.fa.fai $2.percent.bw
#bedGraphToBigWig $3.percent.bedGraph /project/rhie_130/suhn/shared/genomes_new/$5/sequence.fa.fai $3.percent.bw
#rm $2.percent.bedGraph $3.percent.bedGraph

# make samples.txt 
# 2 = hcg 3 = gch
#hcgfile=`realpath $2.percent.bw`
hcgfile=$2
#gchfile=`realpath $3.percent.bw`
gchfile=$3
cd $Dir
ln -sf $hcgfile HCG.bw
ln -sf $gchfile GCH.bw

realpath=`pwd`
hcgfile_ln=${realpath}/HCG.bw
gchfile_ln=${realpath}/GCH.bw

cat << EOF > $Dir/Samples.txt 
$hcgfile_ln	.	percentage
$gchfile_ln	.	percentage
EOF
# make twosamples XXX why do we need this? we shouldn't
cat << EOF > $Dir/Twosamples.txt
$hcgfile_ln	.	percentage
$gchfile_ln	.	percentage
$hcgfile_ln	.	percentage
$gchfile_ln	.	percentage
EOF
# note: lengends and prefixs are not typos!
echo $Dir
echo $Prefix
echo $L1
echo $N1
echo $Name
# density plot #
perl /project/rhie_130/suhn/zexunwu/Bistools_script/alignWigToBed.SR.pl --density_bar --enrich_max 4.0 --result_dir $Dir --prefixs $Prefix1 --locs $L1 --category_names $N1 --sample_names $Name --experiment_names Methylation --experiment_names Accessibility --rep_num_experiments 1 --rep_num_experiments 1 Samples.txt
# average plots #
perl /project/rhie_130/suhn/zexunwu/Bistools_script/alignWigToBed.SR.pl HCG.bw GCH.bw --locs $L1 --average --prefixs $Prefix1 --plot_x_axis_scale 1000 --data_matrix_scale 1200 --bin_size 20 --bin_size_align 1 --plot_x_axis_scale 1000 --result_dir $Dir --smooth --colors black --colors green --lengends HCG --lengends GCH
# heatmap
perl /project/rhie_130/suhn/zexunwu/Bistools_script/alignWigToBed.DKO_paper_version.ranjani.SR.20220228.pl --prefixs $Prefix1 --heatmap_with_reps --data_matrix_scale 1200 --bin_size 20 --bin_size_align 1 --plot_x_axis_scale 1000 --result_dir $Dir --locs $L1 --category_names $N1 --experiment_names Methylation  --experiment_names Accessibility --experiment_names Methylation2  --experiment_names Accessibility2  --heatmap_col blue2yellow --heatmap_col white2darkgreen --heatmap_col blue2yellow --heatmap_col white2darkgreen --heatmap_regionToCluster_low -500 --heatmap_regionToCluster_high 500 --heatmap_cluster 2 --multiSampleClustering 4 Twosamples.txt

