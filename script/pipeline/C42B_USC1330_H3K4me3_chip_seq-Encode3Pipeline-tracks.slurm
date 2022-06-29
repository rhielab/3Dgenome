#!/bin/bash -e
#SBATCH --job-name=C42B_USC1330_H3K4me3_chip_seq-Encode3Pipeline-tracks
#SBATCH --ntasks=1 --cpus-per-task=2 --mem-per-cpu=64GB
#SBATCH --time=47:59:59
#SBATCH --output=C42B_USC1330_H3K4me3_chip_seq-Encode3Pipeline-tracks.out
#SBATCH --error=C42B_USC1330_H3K4me3_chip_seq-Encode3Pipeline-tracks.err
#

module purge
module load gcc/8.3.0
module load anaconda3/2021.05; export PYTHONNOUSERSITE=1
eval "$(conda shell.bash hook)"
conda activate /project/rhie_130/suhn/shared/conda_env/spp
export R_LIBS_USER=/project/rhie_130/suhn/shared/bioinformatic_tools/R4.0_lib
module load igvtools/2.3.98
export PATH="/project/rhie_130/suhn/shared/bioinformatic_tools/bedtools_2.19.1:$PATH"

cd $SLURM_SUBMIT_DIR
cd C42B_USC1330_H3K4me3_chip_seq

FINAL_TA_FILE="C42B_USC1330_H3K4me3_chip_seq.filt.nodup.srt.SE.tagAlign.gz"
PEAK_OUTPUT_DIR=MACS2_PEAKS
CHIP_TA_PREFIX=C42B_USC1330_H3K4me3_chip_seq.filt.nodup.srt.MACS2

# For Fold enrichment signal tracks ============================================ 
# This file is a tab delimited file with 2 columns Col1 (chromosome name), Col2 (chromosome size in bp).
macs2 bdgcmp -t ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}_treat_pileup.bdg -c ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}_control_lambda.bdg --outdir ${PEAK_OUTPUT_DIR} -o ${CHIP_TA_PREFIX}_FE.bdg -m FE
# Remove coordinates outside chromosome sizes (stupid MACS2 bug)
slopBed -i ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}_FE.bdg -g /project/rhie_130/suhn/shared/genomes_new/hg38/hg38.chrom.sizes -b 0 | awk '{if ($3 != -1) print $0}' | /project/rhie_130/suhn/shared/bioinformatic_tools/ucsc/bedClip stdin /project/rhie_130/suhn/shared/genomes_new/hg38/hg38.chrom.sizes ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}.fc.signal.bedgraph
rm -f ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}_FE.bdg
#123456789
# sortBed before bedGraphToBigWig - got error MACS2_PEAKS/SEP044.USC971.H3K27Ac.CON.filt.nodup.srt.MACS2.fc.signal.bedgraph is not case-sensitive sorted at line 45083.  Please use "sort -k1,1 -k2,2n" with LC_COLLATE=C,  or bedSort and try again
# added -k3,3n for overlapping regions
LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}.fc.signal.bedgraph > ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}.fc.signal.sorted.bedgraph

# 123456789
# Convert bedgraph to bigwig using UCSC_Browser_Tools
/project/rhie_130/suhn/shared/bioinformatic_tools/ucsc/bedGraphToBigWig ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}.fc.signal.sorted.bedgraph /project/rhie_130/suhn/shared/genomes_new/hg38/hg38.chrom.sizes ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}.fc.signal.sorted.bw
#rm -f ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}.fc.signal.bedgraph
# Convert bedgraph to tdf file using igv tools
igvtools toTDF -z 7 ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}.fc.signal.sorted.bedgraph ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}.fc.signal.sorted.tdf /project/rhie_130/suhn/shared/genomes_new/hg38/hg38.chrom.sizes 

#=============================== # For -log10(p-value) signal tracks ==============================
# Compute sval = min(no. of reads in ChIP, no. of reads in control) / 1,000,000 
chipReads=$(zcat ${FINAL_TA_FILE} | wc -l | awk '{printf "%f", $1/1000000}');
# 123456789
controlReads=$(zcat /project/rhie_130/suhn/beoungle/encode_chip_seq/c42b/fastq/C42B_control_NA_chip_seq.filt.nodup.srt.PE2SE.tagAlign.gz | wc -l | awk '{printf "%f", $1/1000000}');
sval=$(echo "${chipReads} ${controlReads}" | awk '$1>$2{printf "%f",$2} $1<=$2{printf "%f",$1}');
macs2 bdgcmp -t ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}_treat_pileup.bdg -c ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}_control_lambda.bdg --outdir ${PEAK_OUTPUT_DIR} -o ${CHIP_TA_PREFIX}_ppois.bdg -m ppois -S ${sval}
# Remove coordinates outside chromosome sizes (stupid MACS2 bug)
slopBed -i ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}_ppois.bdg -g /project/rhie_130/suhn/shared/genomes_new/hg38/hg38.chrom.sizes -b 0 | awk '{if ($3 != -1) print $0}' | /project/rhie_130/suhn/shared/bioinformatic_tools/ucsc/bedClip stdin /project/rhie_130/suhn/shared/genomes_new/hg38/hg38.chrom.sizes ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}.pval.signal.bedgraph
rm -rf ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}_ppois.bdg
# 123456789
# added -k3,3n for overlapping regions
LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}.pval.signal.bedgraph > ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}.pval.signal.sorted.bedgraph

# 123456789
# Convert bedgraph to bigwig
/project/rhie_130/suhn/shared/bioinformatic_tools/ucsc/bedGraphToBigWig ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}.pval.signal.sorted.bedgraph /project/rhie_130/suhn/shared/genomes_new/hg38/hg38.chrom.sizes ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}.pval.signal.sorted.bw
#rm -f ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}.pval.signal.bedgraph 
rm -f ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}_treat_pileup.bdg ${peakFile}_control_lambda.bdg
igvtools toTDF -z 7 ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}.pval.signal.sorted.bedgraph ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}.pval.signal.sorted.tdf /project/rhie_130/suhn/shared/genomes_new/hg38/hg38.chrom.sizes
