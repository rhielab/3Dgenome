#!/bin/bash -e
#SBATCH --job-name=C42B_USC1330_H3K4me3_chip_seq-Encode3Pipeline-MACS2
#SBATCH --ntasks=1 --cpus-per-task=4 --mem-per-cpu=32GB
#SBATCH --time=47:59:59
#SBATCH --output=C42B_USC1330_H3K4me3_chip_seq-Encode3Pipeline-MACS2.out
#SBATCH --error=C42B_USC1330_H3K4me3_chip_seq-Encode3Pipeline-MACS2.err
#

module purge
module load gcc/8.3.0
module load anaconda3/2021.05; export PYTHONNOUSERSITE=1
eval "$(conda shell.bash hook)"
conda activate /project/rhie_130/suhn/shared/conda_env/spp
export R_LIBS_USER=/project/rhie_130/suhn/shared/bioinformatic_tools/R4.0_lib
cd $SLURM_SUBMIT_DIR
cd C42B_USC1330_H3K4me3_chip_seq

#=========================
## 2b. Calculate Cross-correlation QC scores

FINAL_TA_FILE="C42B_USC1330_H3K4me3_chip_seq.filt.nodup.srt.SE.tagAlign.gz"

NTHREADS=4

CC_SCORES_FILE="C42B_USC1330_H3K4me3_chip_seq.filt.nodup.sample.15.SE.tagAlign.gz.cc.qc"
CC_PLOT_FILE="C42B_USC1330_H3K4me3_chip_seq.filt.nodup.sample.15.SE.tagAlign.gz.cc.plot.pdf"
# CC_SCORE FILE format
# Filename <tab> numReads <tab> estFragLen <tab> corr_estFragLen <tab> PhantomPeak <tab> corr_phantomPeak <tab> argmin_corr <tab> min_corr <tab> phantomPeakCoef <tab> relPhantomPeakCoef <tab> QualityTag
Rscript --vanilla /project/rhie_130/suhn/shared/bioinformatic_tools/phantompeakqualtools/run_spp_nodups.R -rf -c=C42B_USC1330_H3K4me3_chip_seq.filt.nodup.sample.15.SE.tagAlign.gz -p=${NTHREADS} -filtchr=chrM -savp=${CC_PLOT_FILE} -out=${CC_SCORES_FILE}
#sed -r 's/,[^	]+//g' ${CC_SCORES_FILE} > temp
#mv temp ${CC_SCORES_FILE}
###############################################
#123456789
# without savr and savd - really fast #
read FRAGLEN <<< $(sed -r 's/,[^\t]+//g' $CC_SCORES_FILE | awk 'BEGIN { OFS="\t"} {print $3}')
# sometimes fraglen is nonsense
if [ "$FRAGLEN" = "NA" ] ||  [ $FRAGLEN -lt 50 ]; then
        echo "WARNING: unlikely fraglen $FRAGLEN, using 250"
        FRAGLEN=250
fi
#123456789
PEAK_OUTPUT_DIR=MACS2_PEAKS
mkdir -p $PEAK_OUTPUT_DIR
CHIP_TA_PREFIX="C42B_USC1330_H3K4me3_chip_seq.filt.nodup.srt.MACS2"
##################################################
#========================
# 6. For Histone Marks
# Peak calling and signal tracks using MACSv2 for histone marks

PEAKTYPE=narrow

if [ "$PEAKTYPE" == "narrow" ];
then
    #========= # Generate narrow peaks and preliminary signal tracks ============
    macs2 callpeak -t ${FINAL_TA_FILE} -c /project/rhie_130/suhn/beoungle/encode_chip_seq/c42b/fastq/C42B_control_NA_chip_seq.filt.nodup.srt.PE2SE.tagAlign.gz -f BED -n ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX} -g hs -p 1e-2 --nomodel --shift 0 --extsize ${FRAGLEN} --keep-dup all -B --SPMR
    # Sort by Col8 in descending order and replace long peak names in Column 4 with Peak_<peakRank>
    sort -k 8gr,8gr ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}_peaks.narrowPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' | gzip -c > ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}.narrowPeak.gz
    # remove additional files
    #rm -f ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}_peaks.xls ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}_peaks.narrowPeak ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}_summits.bed
fi 

if [[ "$PEAKTYPE" == "broad" || "$PEAKTYPE" == "gapped" ]]; 
then
    #========= # Generate Broad and Gapped Peaks ================================
    macs2 callpeak -t ${FINAL_TA_FILE} -c /project/rhie_130/suhn/beoungle/encode_chip_seq/c42b/fastq/C42B_control_NA_chip_seq.filt.nodup.srt.PE2SE.tagAlign.gz -n ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX} -g hs -p 1e-2 --broad --nomodel --shift 0 --extsize ${FRAGLEN} --keep-dup all -B --SPMR

    if [ "$PEAKTYPE" == "broad" ];
then
         # Sort by Col8 (for broadPeak) in descending order and replace long peak names in Column 4 with Peak_<peakRank>
         sort -k 8gr,8gr ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}_peaks.broadPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' | gzip -c > ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}.broadPeak.gz
fi

    if [ "$PEAKTYPE" == "gapped" ];
then
         # Sort by Col 14(for gappedPeak) in descending order and replace long peak names in Column 4 with Peak_<peakRank>
         sort -k 14gr,14gr ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}_peaks.gappedPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' | gzip -c > ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}.gappedPeak.gz
fi

    # remove additional files
    #rm -f ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}_peaks.xls ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}_peaks.broadPeak ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}_peaks.gappedPeak ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}_summits.bed
#===========================================
fi
