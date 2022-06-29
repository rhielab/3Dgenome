#!/bin/bash -e
#SBATCH --job-name=C42B_USC1330_H3K4me3_chip_seq-Encode3Pipeline-computeLibraryComplexity-singleEnd
#SBATCH --ntasks=1 --cpus-per-task=2 --mem-per-cpu=64GB
#SBATCH --time=47:59:59
#SBATCH --output=C42B_USC1330_H3K4me3_chip_seq-Encode3Pipeline-computeLibraryComplexity-singleEnd.out
#SBATCH --error=C42B_USC1330_H3K4me3_chip_seq-Encode3Pipeline-computeLibraryComplexity-singleEnd.err
#

module purge
module load gcc/8.3.0
module load python/3.7.6
export PYTHONNOUSERSITE=1
module load openblas/0.3.8 r/4.0.0
export PATH="/project/rhie_130/suhn/shared/bioinformatic_tools/bedtools_2.19.1:$PATH"

cd $SLURM_SUBMIT_DIR
cd C42B_USC1330_H3K4me3_chip_seq

# Compute library complexity
# sort by position and strand

TMPDIR=/scratch/beoungle/chip_tmp_45181735813185/samtools/
FILT_BAM_PREFIX=C42B_USC1330_H3K4me3_chip_seq.filt.srt
FILT_BAM_FILE=${FILT_BAM_PREFIX}.bam

# Obtain unique count statistics
# use bedtools 2.19.1
PBC_FILE_QC="C42B_USC1330_H3K4me3_chip_seq.filt.nodup.srt.pbc.qc"
bedtools bamtobed -i ${FILT_BAM_FILE} | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > ${PBC_FILE_QC}

rm ${FILT_BAM_FILE}
# PBC File output
# TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair
