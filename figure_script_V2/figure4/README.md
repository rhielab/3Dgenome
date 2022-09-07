## Figure4 Regulatory elements and nucleosome-depleted regions (NDRs) that are involved in loops.

This figure is about the regulatory elements (Promoter, Enhancer, Insulator, NDR, repressed regions etc.) involved in loops. We used ChIP-seq data of histone marks (Promoter - H3K4me3, Enhancer - H3K27ac, Repressed - H2K27me3 etc.) as well as CTCF to identify these regions. We processed ChIP-seq data based on the pipeline from ENCODE. Please check https://www.encodeproject.org/chip-seq/histone/
### 4A. Genome browser screenshots of ChIP-seq, NOMe-seq, Hi-C, Micro-C, and RefSeq Genes. 4B. Fractions of regulatory elements that intersect with loop anchors identified from Hi-C 1 billion, Micro-C 1 billion, 2 billion, and 3 billion data. 4C. Numbers of loops belong to different loop categories defined by intersecting the loop anchors with different types of regulatory elements

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
The output is a bunch of loop files with different regulatory elements combination, meaning that the loops with different regulatory elements at two anchors. Additionally, a ```RegSummary.txt``` for the number of loops in different categories will also be created.

### 4D. Comparison of number of promoter-enhancer loops identified from Hi-C 1 billion, Micro-C 1 billion, 2 billion and 3 billion data. 

### 4E. Significance of chromatin interaction (q-value identified by Mustache) for top 5 loop categories. A mean q-value is shown in red. A median q-value is shown in blue.

Here we show an example Rscript ```plot_average_q_value.R``` to generate the q-value violin plot of different chromatin interactions. The input ```.tsv``` files are generated from the ```Loop_Rep_overlapAnalysis.R``` Rscript as described above. Please noted that you need to first generate loop with regulatory element categories, as described above.
