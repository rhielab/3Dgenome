## Figure4 Regulatory elements and nucleosome-depleted regions (NDRs) that are involved in loops.

### (A) Genome browser screenshots of ChIP-seq, NOMe-seq, Hi-C, Micro-C, and RefSeq Genes.  (B) Fractions of regulatory elements that intersect with loop anchors identified from Hi-C 1 billion, Micro-C 1 billion, 2 billion, and 3 billion data. (C) Numbers of loops belong to different loop categories defined by intersecting the loop anchors with different types of regulatory elements

Here we provide an Rscript ```Loop_Rep_overlapAnalysis.R``` to perform the genomic regulatory element distribution analysis for loop data. See below instruction if you are interested to use this script.

First you need to specify the path of regulatory element bed files in script, which include:
- h3k4me3_bed -- path to H3K4me3 or any bed file you used as Promoter.
- h3k27ac_bed -- path to H3K4me3 or any bed file you used as Enhancer.
- CTCF_bed -- path to H3K4me3 or any bed file you used as Insulator.
- ndr_bed -- The bed file you used as NDR.
- h3k27me3_bed -- path to H3K27me3 or any bed file you used as Genomic Repressed region.
- h3k9me3_bed -- path to  H3K9me3 or any bed file you used as HeteroChromatin.

And then you can run the analysis:
```
source('Loop_Rep_overlapAnalysis.R')

My_loop <- read.table(file = '/PATH/TO/LOOPDATA')

loop_reg_analysis(
loop_data = My_loop, 
label = 'OutputPrefix', 
output_path = './output' 
)
```
