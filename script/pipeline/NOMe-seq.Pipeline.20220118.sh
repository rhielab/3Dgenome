#!/bin/bash -e

BIOINFTOOLS=/project/rhie_130/suhn/shared/bioinformatic_tools

if [ $# != 1 ]; then
	echo "usage: $0 filelist"
	echo "Filelist format (no spaces in file names): 1 run per line. Each line contains:"
	echo "GENOME INPUT_NAME OUTPUT_FOLDER INPUT_FASTQS (one if SE, two if PE)"
	exit 1
fi

l=0
cat "$1" | while read GENOME INPUT_NAME OUTPUT_FOLDER INPUT_FASTQS; do
l=$(($l+1))
if [ "$GENOME" = "" ] && [ "$INPUT_NAME" = "" ] && [ "$OUTPUT_FOLDER" = "" ] && [ "$INPUT_FASTQS" = "" ]; then
	continue # empty line
fi
if [ "$GENOME" = "" ] || [ "$INPUT_NAME" = "" ] || [ "$OUTPUT_FOLDER" = "" ] || [ "$INPUT_FASTQS" = "" ]; then
	echo "Missing element in filelist on line $l"
	echo "Filelist format (no spaces in file names): 1 run per line. Each line contains:"
	echo "GENOME INPUT_NAME OUTPUT_FOLDER INPUT_FASTQS (one if SE, two if PE)"
	exit 1
fi

REF="/project/rhie_130/suhn/shared/genomes_new/$GENOME/sequence.fa"
BISCUITASSETS="/project/rhie_130/suhn/shared/genomes_new/$GENOME/biscuit_qc_assets"
BISCUITIDX="/project/rhie_130/suhn/shared/genomes_new/$GENOME/biscuit_index/index"
RGID="$(echo $INPUT_NAME | cut -d _ -f1)"
OUTPUT_PREFIX="$OUTPUT_FOLDER/${GENOME}_$INPUT_NAME"
echo "Genome: $GENOME"
echo "Input name: $INPUT_NAME"
echo "Output folder: $OUTPUT_FOLDER"
echo "Read group ID: $RGID"
echo "Input FASTQ(s): $INPUT_FASTQS"
# check whether to skip the VCF in QC because it would fail
vcfargs="--vcf ${GENOME}_$INPUT_NAME.vcf"
if [ -f "/project/rhie_130/suhn/shared/genomes_new/$GENOME/skip_vcf" ]; then
	echo "Skipping QC functions that require a VCF file because they would fail with this genome"
	vcfargs=
fi
# change some args if SE
SE=0
QCendarg=
if [[ "$INPUT_FASTQS" != *" "* ]]; then
	SE=1
	QCendarg=--single-end
fi
cat << EOF | sbatch
#!/bin/bash -e
#SBATCH --mail-user=$USER@usc.edu --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --partition=oneweek
#SBATCH --time=167:59:59 
#SBATCH --cpus-per-task=8
#SBATCH --mem=59G # 64GB nodes have 59GB max according to CARC
#SBATCH --job-name=$INPUT_NAME
#SBATCH --output=$OUTPUT_PREFIX.log
#SBATCH --error=$OUTPUT_PREFIX.log

module purge
export R_LIBS_USER=$BIOINFTOOLS/R4.0_lib
module load gcc/8.3.0 openblas/0.3.8 r/4.0.3 samtools/1.10 bedtools2/2.27.1 parallel/20190222 picard/2.20.8 htslib/1.10.2
export PATH=$BIOINFTOOLS/biscuit:\$PATH
export PATH=$BIOINFTOOLS/ucsc:\$PATH
mkdir -p "$OUTPUT_FOLDER"
(echo "Read Counts (fastq lines divided by 4)"; echo; 
for i in $INPUT_FASTQS; do printf "%s\t%s\n" "\$i" "\$((\$(zcat -f "\$i" | wc -l)/4))"
done) > $OUTPUT_PREFIX.stats.txt
# No quotes for \$INPUT_FASTQS lets it work with multiple
biscuit align -@ 4 "$BISCUITIDX" $INPUT_FASTQS | samtools sort --write-index -o $OUTPUT_PREFIX.sorted.bam -O BAM -@ 4 
samtools index -@ 8 $OUTPUT_PREFIX.sorted.bam
biscuit pileup -@ 8 -r $REF -N $OUTPUT_PREFIX.sorted.bam -o $OUTPUT_PREFIX.vcf
biscuit vcf2bed -k 1 -t hcg $OUTPUT_PREFIX.vcf > $OUTPUT_PREFIX.HCG.bedGraph
awk -F "\t" 'BEGIN { OFS="\t"} {print \$1,\$2,\$3,\$4}' $OUTPUT_PREFIX.HCG.bedGraph > $OUTPUT_PREFIX.HCG.noCounts.bedGraph
biscuit vcf2bed -k 1 -t gch $OUTPUT_PREFIX.vcf > $OUTPUT_PREFIX.GCH.bedGraph
awk -F "\t" 'BEGIN { OFS="\t"} {print \$1,\$2,\$3,\$4}' $OUTPUT_PREFIX.GCH.bedGraph > $OUTPUT_PREFIX.GCH.noCounts.bedGraph
bedGraphToBigWig $OUTPUT_PREFIX.HCG.noCounts.bedGraph $REF.fai $OUTPUT_PREFIX.HCG.bw
bedGraphToBigWig $OUTPUT_PREFIX.GCH.noCounts.bedGraph $REF.fai $OUTPUT_PREFIX.GCH.bw
gzip -c $OUTPUT_PREFIX.vcf > $OUTPUT_PREFIX.vcf.gz
rm -f $OUTPUT_PREFIX.vcf
biscuit vcf2bed -k 1 -t hcg $OUTPUT_PREFIX.vcf.gz | awk 'BEGIN { OFS="\\t"} {print \$1,\$2,\$3,".",int(\$4*100),\$6=="C"?"+":"-",\$4,\$5}' > $OUTPUT_PREFIX.HCG.6plus2.bed
biscuit vcf2bed -k 1 -t gch $OUTPUT_PREFIX.vcf.gz | awk 'BEGIN { OFS="\\t"} {print \$1,\$2,\$3,".",int(\$4*100),\$6=="C"?"+":"-",\$4,\$5}' > $OUTPUT_PREFIX.GCH.6plus2.bed
awk 'BEGIN { OFS="\t"} {print \$1,\$2,\$3,\$4,\$5,\$6,100*\$7,\$8}' $OUTPUT_PREFIX.GCH.6plus2.bed > $OUTPUT_PREFIX.GCH.6plus2.percent.bed
# make percent-based bigWigs for Bis-plot scripts
awk 'BEGIN { OFS="\t"} {print $1, $2, $3, $4*100}' $OUTPUT_PREFIX.GCH.noCounts.bedGraph > $OUTPUT_PREFIX.GCH.noCounts.percent.bedGraph
awk 'BEGIN { OFS="\t"} {print $1, $2, $3, $4*100}' $OUTPUT_PREFIX.HCG.noCounts.bedGraph > $OUTPUT_PREFIX.HCG.noCounts.percent.bedGraph
bedGraphToBigWig $OUTPUT_PREFIX.GCH.noCounts.percent.bedGraph $REF.fai $OUTPUT_PREFIX.GCH.noCounts.percent.bw
bedGraphToBigWig $OUTPUT_PREFIX.HCG.noCounts.percent.bedGraph $REF.fai $OUTPUT_PREFIX.HCG.noCounts.percent.bw
# convert to TSV to work with aaRon
printf "chr\tposition\t%s.C\t%s.cov\n" "$INPUT_NAME" "$INPUT_NAME" > $OUTPUT_PREFIX.GCH.6plus2.aaRon.tsv
awk 'BEGIN { OFS="\t"} {print \$1, \$2, \$8*\$7/100, \$8}' $OUTPUT_PREFIX.GCH.6plus2.bed >> $OUTPUT_PREFIX.GCH.6plus2.aaRon.tsv
Rscript $BIOINFTOOLS/aaRon/run_aaRon.R "$GENOME" "$OUTPUT_PREFIX" "$INPUT_NAME" "$OUTPUT_PREFIX.GCH.6plus2.aaRon.tsv" 
# QC starts here
picard MarkDuplicates CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT METRICS_FILE=$OUTPUT_PREFIX.sorted.mdups.metrics.txt READ_NAME_REGEX=null INPUT=$OUTPUT_PREFIX.sorted.bam OUTPUT=$OUTPUT_PREFIX.sorted.mdups.bam
samtools flagstat -@ 7 $OUTPUT_PREFIX.sorted.mdups.bam > $OUTPUT_PREFIX.sorted.mdups.flagstat.txt # in flagstat, for some reason, -@ is the number of ADDITIONAL threads
# Obtain unique count statistics
bedtools bamtobed -i $OUTPUT_PREFIX.sorted.mdups.bam | awk 'BEGIN{OFS="\\t"}{print \$1,\$2,\$3,\$6}' | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} (\$1==1){m1=m1+1} (\$1==2){m2=m2+1} {m0=m0+1} {mt=mt+\$1} END{printf "%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > $OUTPUT_PREFIX.sorted.mdups.pbc_qc.txt
# PBC File output:
# TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair
TOTAL_READ_COUNT="\$(awk 'BEGIN{a=0}(NR>2){a+=\$2}END{print a}' $OUTPUT_PREFIX.stats.txt)"
if [ $SE = 1 ]; then
	UNIQ_COUNT="\$(cat $OUTPUT_PREFIX.sorted.mdups.flagstat.txt | awk 'NR==1{map=\$1}END{print map}')"
else
	UNIQ_COUNT=\$(cat $OUTPUT_PREFIX.sorted.mdups.flagstat.txt | grep 'read1' | cut -d ' ' -f 1)
fi
PBC_STATS="\$(awk '{print \$1" "\$2" "\$5" "\$6" "\$7}' $OUTPUT_PREFIX.sorted.mdups.pbc_qc.txt)"
echo >> $OUTPUT_PREFIX.stats.txt
echo -e "TOTAL_READS (from fastq)\tUNIQ_READS\tTOTAL_READS(PBC_STATS)\tDISTINCT_READS(PBC_STATS)\tNRF=DISTINCT/TOTAL(PBC_STATS)\tPBC1=OnePair/Distinct\tPBC2=OnePair/TwoPair" >> $OUTPUT_PREFIX.stats.txt
echo -e "\$TOTAL_READ_COUNT\t\$UNIQ_COUNT\t\$PBC_STATS" >> $OUTPUT_PREFIX.stats.txt
#XXXX 
#picard CollectAlignmentSummaryMetrics INPUT=$OUTPUT_PREFIX.sorted.mdups.bam OUTPUT=$OUTPUT_PREFIX.sorted.mdups.CollectAlignmentSummaryMetrics.txt IS_BISULFITE_SEQUENCED=true VALIDATION_STRINGENCY=SILENT REFERENCE_SEQUENCE=$REF
echo "CollectAlignmentSummaryMetrics doesn't work with Biscuit right now, skipping" | tee $OUTPUT_PREFIX.sorted.mdups.CollectAlignmentSummaryMetrics.txt
picard AddOrReplaceReadGroups CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=1000000 INPUT=$OUTPUT_PREFIX.sorted.mdups.bam OUTPUT=$OUTPUT_PREFIX.sorted.mdups.w.rg.bam RGID=$RGID RGLB=$RGID RGPL="Illumina" RGPU=$RGID RGSM=$RGID RGCN="Rhie Lab" RGDS="From file $(basename $OUTPUT_PREFIX.sorted.mdups.w.rg.bam) on $(date '+%Y-%m-%d')"
picard CollectInsertSizeMetrics I=$OUTPUT_PREFIX.sorted.mdups.w.rg.bam O=$OUTPUT_PREFIX.sorted.mdups.w.rg.CollectInsertSizeMetrics.txt H=$OUTPUT_PREFIX.sorted.mdups.w.rg.CollectInsertSizeMetrics.pdf VALIDATION_STRINGENCY=SILENT
# check coverage projection metrics
$BIOINFTOOLS/preseq/preseq lc_extrap -v -s 10000000 -e 3000000000 -bam $OUTPUT_PREFIX.sorted.mdups.w.rg.bam > $OUTPUT_PREFIX.sorted.mdups.w.rg.CoverageProjection.metric.txt || (echo "preseq lc_extrap failed - too few reads?" | tee $OUTPUT_PREFIX.sorted.mdups.w.rg.CoverageProjection.metric.txt)
# run Biscuit QC. It segfaults if the paths are too long, so cd into the output directory
cd "$OUTPUT_FOLDER"
$BIOINFTOOLS/biscuit/scripts/QC.sh $QCendargs $vcfargs --outdir ${GENOME}_$INPUT_NAME.biscuit.qc $BISCUITASSETS $REF $RGID ${GENOME}_$INPUT_NAME.sorted.mdups.w.rg.bam
# if totalReadConversionRate fails, try the old way
# based on https://github.com/genome/docker-biscuit/blob/46231ab60727494bce7aeb5fdbf125361daf908e/Bisulfite_QC_bisulfiteconversion.sh
# samtools view command updated from https://github.com/huishenlab/biscuit/blob/013e6b0e839e231d8e7f7c08984da6adf21eee25/scripts/QC.sh#L139
if fgrep -q -- "nan" "${GENOME}_$INPUT_NAME.biscuit.qc/${RGID}_totalReadConversionRate.txt"; then
	echo "totalReadConversionRate.txt generation failed, trying old way"
	echo -e "BISCUITqc Conversion Rate by Read Average Table" > "${GENOME}_$INPUT_NAME.biscuit.qc/${RGID}_totalReadConversionRate.txt"
	samtools view -F 0x500 -q 40 -h ${GENOME}_$INPUT_NAME.sorted.mdups.w.rg.bam | biscuit bsconv -p $REF - | awk '{for(i=1;i<=8;++i) a[i]+=$i;}END{print "CpA\tCpC\tCpG\tCpT"; print a[1]/(a[1]+a[2])"\t"a[3]/(a[3]+a[4])"\t"a[5]/(a[5]+a[6])"\t"a[7]/(a[7]+a[8]);}' >> "${GENOME}_$INPUT_NAME.biscuit.qc/${RGID}_totalReadConversionRate.txt"
fi
# additional QC commands from Biscuit site
biscuit bsconv $REF ${GENOME}_$INPUT_NAME.sorted.mdups.w.rg.bam ${GENOME}_$INPUT_NAME.sorted.mdups.w.rg.bsconv_tags.bam
biscuit bsstrand $REF ${GENOME}_$INPUT_NAME.sorted.mdups.w.rg.bam > ${GENOME}_$INPUT_NAME.biscuit.qc/${RGID}_strand_table.txt # replace original one because it has more detail
biscuit cinread -t hcg -o ${GENOME}_$INPUT_NAME.biscuit.qc/${RGID}_cytosine-read_pairs_HCG.txt $REF ${GENOME}_$INPUT_NAME.sorted.mdups.w.rg.bam
biscuit cinread -t gch -o ${GENOME}_$INPUT_NAME.biscuit.qc/${RGID}_cytosine-read_pairs_GCH.txt $REF ${GENOME}_$INPUT_NAME.sorted.mdups.w.rg.bam
# generate read depth histogram PDFs
Rscript --vanilla /project/rhie_130/suhn/shared/NOMe-seq.Pipeline/ReadDepthHistogram.R $OUTPUT_PREFIX.HCG.6plus2.bed
Rscript --vanilla /project/rhie_130/suhn/shared/NOMe-seq.Pipeline/ReadDepthHistogram.R $OUTPUT_PREFIX.GCH.6plus2.bed
EOF
echo # to separate runs
done
