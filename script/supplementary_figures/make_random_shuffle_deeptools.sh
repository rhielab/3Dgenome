module load gcc/8.3.0
module load anaconda3

eval "$(conda shell.bash hook)"

conda activate deeptools

TMPDIR=/project/rhie_130/suhn/beoungle/combined/temp

for i in {1..10};
do
prefixing=$i
second="H3K27ac_peak_in_Micro_C_3_Bil.chrndr.bed"
real_prefix="${second/chrndr/$prefixing}"
shuf -n 1000 H3K27ac_peak_in_Micro_C_3_Bil.chrndr.bed | awk -v OFS='\t' '{print$1,$2,$3}' > $real_prefix 
computeMatrix reference-point --referencePoint center  -b 1000 -a 1000 -S /project/rhie_130/suhn/beoungle/combined/mapped_C42B_Novagene_all.PT.bigwig -R $real_prefix -p 16 -o $real_prefix.1kb.RPKM.gz
plotHeatmap -m $real_prefix.1kb.RPKM.gz -o $real_prefix.pdf --legendLocation none
done

for i in {1..10};
do
prefixing=$i
second="NDR_no_feature_in_Micro_C_3_Bil.chr.bed"
real_prefix="${second/chr/$prefixing}"
shuf -n 10000 NDR_no_feature_in_Micro_C_3_Bil.chr.bed | awk -v OFS='\t' '{print$1,$2,$3}' > $real_prefix
computeMatrix reference-point --referencePoint center  -b 1000 -a 1000 -S /project/rhie_130/suhn/beoungle/combined/mapped_C42B_Novagene_all.PT.bigwig -R $real_prefix -p 16 -o $real_prefix.1kb.RPKM.gz
plotHeatmap -m $real_prefix.1kb.RPKM.gz -o $real_prefix.pdf --legendLocation none
done

for i in {1..10};
do
prefixing=$i
second="CTCF_motif_peak_in_Micro_C_3_Bil.chr.bed"
real_prefix="${second/chr/$prefixing}"
shuf -n 30000 CTCF_motif_peak_in_Micro_C_3_Bil.chr.bed | awk -v OFS='\t' '{print$1,$2,$3}' > $real_prefix
computeMatrix reference-point --referencePoint center  -b 1000 -a 1000 -S /project/rhie_130/suhn/beoungle/combined/mapped_C42B_Novagene_all.PT.bigwig -R $real_prefix -p 16 -o $real_prefix.1kb.RPKM.gz
plotHeatmap -m $real_prefix.1kb.RPKM.gz -o $real_prefix.pdf --legendLocation none
done

for i in {1..10};
do
prefixing=$i
second="TSS_in_Micro_C_3_Bil.chr.bed"
real_prefix="${second/chr/$prefixing}"
shuf -n 10000 TSS_in_Micro_C_3_Bil.chr.bed | awk -v OFS='\t' '{print$1,$2,$3}' > $real_prefix
computeMatrix reference-point --referencePoint center  -b 1000 -a 1000 -S /project/rhie_130/suhn/beoungle/combined/mapped_C42B_Novagene_all.PT.bigwig -R $real_prefix -p 16 -o $real_prefix.1kb.RPKM.gz
plotHeatmap -m $real_prefix.1kb.RPKM.gz -o $real_prefix.pdf --legendLocation none
done


