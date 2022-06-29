#!/bin/bash -e
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=32GB
#SBATCH --time=23:59:59
#SBATCH --mail-user=zexunwu@usc.edu
#SBATCH --mail-type=all
#SBATCH --output=chicago.out
#SBATCH --error=chicago.err
#SBATCH --partition=rhie

#Bam file you got from matrices at the end, with PT name.
#You can obtain bait, probe, and other information from this page as well https://github.com/dovetail-genomics/capture/blob/main/docs/source/data_sets.rst
#This script uses everything from https://github.com/dovetail-genomics/capture/tree/main/docs/source

bam=/project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/Chicago/input/PC-UNI1945.bam
outputpath=/project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/Chicago/output2
baits=/project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/Chicago/input/baits_v1.0.bed.bgz

fname=$(basename $bam)
id=${fname%.bam}

res=5

cd $outputpath
#Clean up the bam file so it is compatbile with Chicago later on.
module load samtools
module load bedtools2
module load gcc
module load r

#remove -n in samtools sort because the name-sorted bam files are not compatiable with samtools index. - Tim
#samtools view -@ 8 -Shu -F 2048 ${bam} | \
#    samtools sort -n -@ 8 -o ${id}-cleanedUp.bam

#samtools index -@ 8 ${id}-cleanedUp.bam - comment this out because samtools index doesn't work with name-sorted bam file

samtools view -H ${id}-cleanedUp.bam | \
    awk -v OFS='\t' '/^@SQ/ && !($2 ~ /:(chr|"")M/) {split($2,chr,":");split($3,ln,":");print chr[2],ln[2]}' | \
    sort -V -k1,1 \
	 > chr_size.tsv


## prep4Chicago.R  ${baits} $(( ${res} * 1000 )) ${id}-cleanedUp.bam

#Creating baitmap that matches resolution that we are putting in.
#Why use below command to generate files?
zcat $baits \
    | awk 'OFS="\t" {$1="chr"$1; print}' \
    | bedtools sort -i stdin \
    | cut -f1,2,3 \
    | bedtools slop -g chr_size.tsv -b 120 -i stdin \
    | awk -F'\t' 'NR>0{$0=$0"\t""bait_"NR} 1' \
    | bedtools merge -i stdin -c 4 -o collapse -delim "_" \
	       > human_pooled_baits120bp.bed

bedtools makewindows -g chr_size.tsv -w $(( ${res} * 1000 )) > genome.${res}kb.bed

bedtools subtract -a genome.${res}kb.bed -b human_pooled_baits120bp.bed > ${res}kb_sub_probe.bed


cat human_pooled_baits120bp.bed \
    <(awk '{print $1"\t"$2"\t"$3"\t""label"}' ${res}kb_sub_probe.bed) \
    | bedtools sort -i stdin \
    | awk -F'\t' 'NR>0{$0=$0"\t"NR} 1' \
	  > ${res}kb_OEnBaits.bed

awk '{print $1"\t"$2"\t"$3"\t"$5}' ${res}kb_OEnBaits.bed  > ${res}kb.rmap

awk '{if ($4 != "label") print $1"\t"$2"\t"$3"\t"$5"\t"$4}' ${res}kb_OEnBaits.bed  > ${res}kb.baitmap
    



#It's design files; typo on Desing from dovetail. Dovetail also provide design files, but not for higher resolution at 1kb or 2kb, so we need to generate them anyway.
python ../makeDesignFiles_py3.py \
       --minFragLen 75 \
       --maxFragLen 7000 \
       --maxLBrownEst 1000000 \
       --binsize 20000 \
       --rmapfile ${res}kb.rmap \
       --baitmapfile ${res}kb.baitmap \
       --outfilePrefix ${res}kDesingFiles

#Refer to this section https://github.com/dovetail-genomics/capture/blob/main/docs/source/interactions.rst
#Convert bam file into chicago compatible file with baitmap and rmap information.
../bam2chicago.sh ${id}-cleanedUp.bam  ${res}kb.baitmap  ${res}kb.rmap  chinput_${id}

#Run chicago to identify chromatin interactions. Make sure that ID name matches output name here, or they won't run properly. (i.e. PC-UNI146.bam, needs PC-UNI1946 at the end, and chinput_PC-UNI146/chinput_PC-UNI1946.chinput for result section. You can change design-dir papamter to where design files will actually be at.
Rscript ../runChicago.R \
	--design-dir . \
	--cutoff 5 \
	--export-format interBed,washU_text,seqMonk,washU_track \
	chinput_${id}/chinput_${id}.chinput \
	${id}
