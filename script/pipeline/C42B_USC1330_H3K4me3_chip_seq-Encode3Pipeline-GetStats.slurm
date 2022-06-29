#!/bin/bash -e
#SBATCH --job-name=C42B_USC1330_H3K4me3_chip_seq-Encode3Pipeline-GetStats
#SBATCH --ntasks=1 --cpus-per-task=2 --mem-per-cpu=64GB
#SBATCH --time=47:59:59
#SBATCH --output=C42B_USC1330_H3K4me3_chip_seq-Encode3Pipeline-GetStats.out
#SBATCH --error=C42B_USC1330_H3K4me3_chip_seq-Encode3Pipeline-GetStats.err
#

module purge

cd $SLURM_SUBMIT_DIR
cd C42B_USC1330_H3K4me3_chip_seq

TOTAL_READ_COUNT="$(zcat -f /project/rhie_130/suhn/beoungle/encode_chip_seq/c42b/fastq/C42B_H3K4me3_Rep1_USC1330_K4me3_C42B_S6_L001_R1_001.fastq.gz | echo $((`wc -l`/4)))"
UNIQ_COUNT="$(cat C42B_USC1330_H3K4me3_chip_seq.filt.nodup.srt.flagstat.qc | awk 'NR==1{map=$1}END{print map}')"
PBC_STATS="$(cat C42B_USC1330_H3K4me3_chip_seq.filt.nodup.srt.pbc.qc | awk '{print $1" "$2" "$5" "$6" "$7}')"
CCQC_STATS="$(cat *tagAlign.gz.cc.qc  | awk '{print $3" "$9" "$10}' | awk '{split($1,a,","); print a[1]" "$2" "$3}')" 
echo "TOTAL_READS (from fastq)	UNIQ_READS	TOTAL_READS(PBC_STATS)	DISTINCT_READS(PBC_STATS)	NRF=DISTINCT/TOTAL(PBC_STATS)	PBC1=OnePair/Distinct	PBC2=OnePair/TwoPair	EST_FRAG_LENGTH	PHANTOM_PEAK_COEF(CCQC_STATS)	REL_PHANTOM_PEAK_COEF(CCQC_STATS)" > C42B_USC1330_H3K4me3_chip_seq-PipelineStats
echo "$TOTAL_READ_COUNT	$UNIQ_COUNT	$PBC_STATS	$CCQC_STATS" >> C42B_USC1330_H3K4me3_chip_seq-PipelineStats

rm -rf "/scratch/beoungle/chip_tmp_45181735813185"

