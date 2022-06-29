#!/bin/bash -e
#SBATCH --job-name=C42B_USC1330_H3K4me3_chip_seq-Encode3Pipeline-bwa-singleEnd
#SBATCH --ntasks=1 --cpus-per-task=4 --mem-per-cpu=32GB
#SBATCH --time=47:59:59
#SBATCH --output=C42B_USC1330_H3K4me3_chip_seq-Encode3Pipeline-bwa-singleEnd.out
#SBATCH --error=C42B_USC1330_H3K4me3_chip_seq-Encode3Pipeline-bwa-singleEnd.err
#
cd $SLURM_SUBMIT_DIR
cd C42B_USC1330_H3K4me3_chip_seq

module purge
module load gcc/8.3.0
module load python/3.7.6
export PYTHONNOUSERSITE=1
module load openblas/0.3.8 r/4.0.0
export PATH="/project/rhie_130/suhn/shared/bioinformatic_tools/bedtools_2.19.1:$PATH"
module load samtools/1.10
module load bwa/0.7.17
module load picard/2.20.8


mkdir -p /scratch/beoungle/chip_tmp_45181735813185/samtools;
mkdir -p /scratch/beoungle/chip_tmp_45181735813185/picard;

#===========================
## 1a. Read alignment (BWA)
###
# Run BWA alignment
###

# works fine both with gz file and unzipped file
FASTQ_FILE_1=/project/rhie_130/suhn/beoungle/encode_chip_seq/c42b/fastq/C42B_H3K4me3_Rep1_USC1330_K4me3_C42B_S6_L001_R1_001.fastq.gz

NTHREADS=4

# map reads to create raw SAM file 
SAI_FILE_1="C42B_USC1330_H3K4me3_chip_seq.sai"
RAW_BAM_PREFIX="C42B_USC1330_H3K4me3_chip_seq.raw.srt"
RAW_BAM_FILE="${RAW_BAM_PREFIX}.bam" 
RAW_BAM_FILE_MAPSTATS="${RAW_BAM_PREFIX}.flagstat.qc" # QC File

bwa aln -q 5 -l 32 -k 2 -t ${NTHREADS} /project/rhie_130/suhn/shared/genomes_new/hg38/bwa_index/index ${FASTQ_FILE_1} > ${SAI_FILE_1}
bwa samse /project/rhie_130/suhn/shared/genomes_new/hg38/bwa_index/index ${SAI_FILE_1} ${FASTQ_FILE_1} | samtools view -Su - | samtools sort -o ${RAW_BAM_FILE}
rm ${SAI_FILE_1}
samtools flagstat ${RAW_BAM_FILE} > ${RAW_BAM_FILE_MAPSTATS}

#===============================
## 1b. Post-alignment filtering 
###
# Run Post-alignment filtering
###

FILT_BAM_PREFIX="C42B_USC1330_H3K4me3_chip_seq.filt.srt"
FILT_BAM_FILE="${FILT_BAM_PREFIX}.bam"
MAPQ_THRESH=30

samtools view -F 1804 -q ${MAPQ_THRESH} -b ${RAW_BAM_FILE} > ${FILT_BAM_FILE}

# Mark duplicates 
TMP_FILT_BAM_FILE="${FILT_BAM_PREFIX}.dupmark.bam"

DUP_FILE_QC="${FILT_BAM_PREFIX}.dup.qc" # QC file
picard MarkDuplicates INPUT=${FILT_BAM_FILE} OUTPUT=${TMP_FILT_BAM_FILE} METRICS_FILE=${DUP_FILE_QC} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false TMP_DIR=/scratch/beoungle/chip_tmp_45181735813185/picard
mv ${TMP_FILT_BAM_FILE} ${FILT_BAM_FILE}

# remove duplicates and index final position sorted BAM
FINAL_BAM_FILE="C42B_USC1330_H3K4me3_chip_seq.filt.nodup.srt.bam" # To be stored 
FINAL_BAM_INDEX_FILE="C42B_USC1330_H3K4me3_chip_seq.filt.nodup.srt.bai" # To be stored 
FINAL_BAM_FILE_MAPSTATS="C42B_USC1330_H3K4me3_chip_seq.filt.nodup.srt.flagstat.qc" # QC file
samtools view -F 1804 -b ${FILT_BAM_FILE} > ${FINAL_BAM_FILE}

# Index Final BAM file
samtools index ${FINAL_BAM_FILE} ${FINAL_BAM_INDEX_FILE}
samtools flagstat ${FINAL_BAM_FILE} > ${FINAL_BAM_FILE_MAPSTATS}
