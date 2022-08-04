#!/bin/bash
#1st input - path to loop data
#2nd input - output path

file_name=$(basename $1)
OUTPUTPATH="${2:-$PWD}"
prefixing=$OUTPUTPATH/$file_name
summary=$OUTPUTPATH/loopDistance_summary.txt

N_total=$(cat $1 | wc -l)
echo '#Summary of loops with different distance category. The distance is calculated using 5th column minus 2nd column' > $summary
echo '#InputFiles:'$(realpath $1) >> $summary
echo -e 'Category\tNumber\tRatio' >> $summary
echo -e 'Total\t'$N_total'\t'1 >> $summary

#For loop <= 200kb
second="less_200kb.bedpe"
real_prefix="${prefixing/tsv/$second}"
awk '{if ($5-$2 <= 200000) print}' $1 > $real_prefix
length=$(cat $real_prefix | wc -l)
ratio=$(echo "scale=6; $length/$N_total" | bc | awk '{printf "%f", $0}')
echo -e 'less_200kb\t'$length'\t'$ratio >> $summary

#For loop 200kb - 400kb
second="200kb_400kb.bedpe"
real_prefix="${prefixing/tsv/$second}"
awk '{if ( $5-$2 > 200000 && $5-$2 <= 400000 ) print}' $1 > $real_prefix
length=$(cat $real_prefix | wc -l)
ratio=$(echo "scale=6; $length/$N_total" | bc | awk '{printf "%f", $0}')
echo -e '200kb_400kb\t'$length'\t'$ratio >> $summary

#For loop 400kb - 600kb
second="400kb_600kb.bedpe"
real_prefix="${prefixing/tsv/$second}"
awk '{if ( $5-$2 > 400000 && $5-$2 <= 600000 ) print}' $1 > $real_prefix
length=$(cat $real_prefix | wc -l)
ratio=$(echo "scale=6; $length/$N_total" | bc | awk '{printf "%f", $0}')
echo -e '400kb_600kb\t'$length'\t'$ratio >> $summary

#For loop 600kb - 800kb
second="600kb_800kb.bedpe"
real_prefix="${prefixing/tsv/$second}"
awk '{if ( $5-$2 > 600000 && $5-$2 <= 800000 ) print}' $1 > $real_prefix
length=$(cat $real_prefix | wc -l)
ratio=$(echo "scale=6; $length/$N_total" | bc | awk '{printf "%f", $0}')
echo -e '600kb_800kb\t'$length'\t'$ratio >> $summary

#For loop 800kb - 1Mb
second="800kb_1Mb.bedpe"
real_prefix="${prefixing/tsv/$second}"
awk '{if ( $5-$2 > 800000 && $5-$2 <= 1000000 ) print}' $1 > $real_prefix
length=$(cat $real_prefix | wc -l)
ratio=$(echo "scale=6; $length/$N_total" | bc | awk '{printf "%f", $0}')
echo -e '800kb_1Mb\t'$length'\t'$ratio >> $summary

#For loop larger than 1Mb
second="greater_1MB.bedpe"
real_prefix="${prefixing/tsv/$second}"
awk '{if ( $5-$2 > 1000000) print}' $1 > $real_prefix
length=$(cat $real_prefix | wc -l)
ratio=$(echo "scale=6; $length/$N_total" | bc | awk '{printf "%f", $0}')
echo -e 'greater_1Mb\t'$length'\t'$ratio >> $summary
