#!/bin/bash -e
#SBATCH --job-name=C42B_USC1330_H3K4me3_chip_seq-Encode3Pipeline-tagAlign-singleEnd
#SBATCH --ntasks=1 --cpus-per-task=2 --mem-per-cpu=64GB
#SBATCH --time=47:59:59
#SBATCH --output=C42B_USC1330_H3K4me3_chip_seq-Encode3Pipeline-tagAlign-singleEnd.out
#SBATCH --error=C42B_USC1330_H3K4me3_chip_seq-Encode3Pipeline-tagAlign-singleEnd.err
#

module purge
module load gcc/8.3.0
module load python/3.7.6
export PYTHONNOUSERSITE=1
module load openblas/0.3.8 r/4.0.0
export PATH="/project/rhie_130/suhn/shared/bioinformatic_tools/bedtools_2.19.1:$PATH"

cd $SLURM_SUBMIT_DIR
cd C42B_USC1330_H3K4me3_chip_seq

#===========================
### 2a. convert BAM to tagAlign (BED 3+3 format)
##

# bedtools 2.19.1, gawk shuf
# module add bedtools/2.19.1

# Create tagAlign file 
FINAL_BAM_FILE="C42B_USC1330_H3K4me3_chip_seq.filt.nodup.srt.bam"
FINAL_TA_FILE="C42B_USC1330_H3K4me3_chip_seq.filt.nodup.srt.SE.tagAlign.gz"

bedtools bamtobed -i ${FINAL_BAM_FILE} | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' | gzip -c > ${FINAL_TA_FILE}

# Subsample tagAlign file
NREADS=15000000
#SUBSAMPLED_TA_FILE="C42B_USC1330_H3K4me3_chip_seq.filt.nodup.sample.$(($NREADS/1000000)).SE.tagAlign.gz"
zcat ${FINAL_TA_FILE} | grep -v "chrM" | shuf -n ${NREADS} | gzip -c > C42B_USC1330_H3K4me3_chip_seq.filt.nodup.sample.15.SE.tagAlign.gz
