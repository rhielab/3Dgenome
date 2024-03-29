## *__Rhie Lab - 3D genome project__*

## Citation
Lee, B.H., Wu, Z. & Rhie, S.K. Characterizing chromatin interactions of regulatory elements and nucleosome positions, using Hi-C, Micro-C, and promoter capture Micro-C. Epigenetics & Chromatin 2022 Dec 21;15(1):41. https://doi.org/10.1186/s13072-022-00473-4 PMID: 36544209

## Datasets
All of Hi-C, Micro-C, Capture Micro-C, ChIP-seq, NOMe-seq, and RNA-seq datasets used for the above study can be found from the GEO GSE205000 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205000

## Analysis
Here we describe the analysis method we used to generate main figures of the above manuscript (PMID: 36544209). 

All of the in-house scripts we used can be found from [/script_new](https://github.com/rhielab/3Dgenome/tree/main/script_new). 

If you use the datasets and scripts we generated for your study, *__please cite (PMID: 36544209)__* to your manuscript.

We processed Hi-C, Micro-C, and Capture Micro-C data using the 4DN Hi-C data processing pipeline (https://data.4dnucleome.org/resources/data-analysis/hi_c-processing-pipeline) to get all of the necessary files .hic, .bam, .cool, .mcool required for different downstream softwares.

We processed ChIP-seq data using the ENCODE ChIP-seq pipeline (https://www.encodeproject.org/chip-seq/histone/) to annnotate regulatory elements (Promoter, Enhancer, Insulator, NDR, Repressed Regions etc.) 

We processed NOMe-seq data using BSMAP (Xi and Li, BMC Bioinformatics, 2009 PMID:19635165) to align .fastq file to bisulfile-converted genome and Bis-SNP (Liu et al, Genome Biol, 2012 PMID:22784381) to identify meethylation status of CpG and GpC sites from BAM file.

## Figure 1. Comparison of Hi-C and Micro-C data
### 1A. Experimental methods of Hi-C and Micro-C 

### 1B. Chromatin interaction heatmaps of Hi-C and Micro-C data
We used cooltools (https://github.com/open2c/cooltools) to generate the interaction heatmap. We used a Python script ```cooltools_heatmap.py``` to run cooltools. The required input file to generate this plot is ```.mcool``` file. Below is the example code:

To run the script, please make sure that you have installed cooltools correctly. 
- 1st - mcool file
- 2nd - output path
- 3rd - prefix for the output
```
python cooltools_heatmap.py \
input/EXAMPLE.mcool \
output/ \
output_prefix
```

### 1C. Venn diagram of TADs identified from Hi-C and Micro-C

We used an Rscript ```topdom_comparison.R``` to perform the overlap analysis using the TAD data generated from TopDom (Shin et al, Nucleic Acids Res, 2016 PMID:26704975). The inputs are the TAD bed files generated by TopDom and the script will output the shared TAD number. Below is example usage:

Load the function. Make sure you have libraries ```dplyr``` and ```fuzzyjoin``` available in your R environment.

```
source('topdom_comparison.R')
```

Then load the TopDom derived bed files:
```
data1 <- read.table('topdom_1.bed')
data2 <- read.table('topdom_2.bed')
```

Run the analysis. Specify the resolution and labels as needed. Deafult for labels are "topdom_bed1" and "topdom_bed2".
```
final <- topdom_comparison(topdom_bed1 = data1, topdom_bed2 = data2, 
 label_1='C42B.1B', label_2 = 'C42B.HiC', resolution = 50000)
```

The output is a list containing four tables. The first table is the overlapped TADs in bed1; second table is the overlapped TADs in bed2; third table is the bed1 with comparison 0/1 matrix; fourth table is bed2 with comparison 0/1 matrix. To check:
```
head(final)
```

For detailed instructions for TopDom, please refer to GitHub page: https://github.com/HenrikBengtsson/TopDom.

### 1D. Triangular heatmaps of Hi-C and Micro-C near chr1p32 region

The triangular heatmaps are generated using ```.hic``` files in [UCSC genome browser](https://genome.ucsc.edu). Please refer to 1C for TAD calling.

### 1E. Average chromatin interaction signals at shared loops and unique loops

We used a Python script ```cooltools_pileup.py``` to generate pileup plot using coolstools. The required input file format is ```.mcool```. Below is the example code:

To run ```cooltools_pileup.py```, cooltools conda environment is required
- 1st - mcool file
- 2nd - loop file 
- 3rd - output path
- 4th - output prefix
```
python cooltools_pileup.py \
EXAMPLE.mcool \
EXAMPLE-loop.bedpe \
./output \
output_prefix
```

### 1F. Triangular heatmaps of Hi-C and Micro-C near chr7p14 region

Please refer to 1D for heatmaps generation. Loops are called by Mustache (Ardakany et al, Genome Biol, 2020 PMID:32998764), which takes ```.hic``` file as input.  Please refer to GitHub page of Mustache for more detailed instructions: https://github.com/ay-lab/mustache.

## Figure 2. Comparison of Micro-C data in different read depth sequencing

### 2A. Chromatin interaction heatmaps of Micro-C 1 billion, 2 billion, and 3 billion data

The plots are generated by using cooltools. 

### 2B. Numbers of loops identified at different resolutions from Hi-C 1 billion, Micro-C 1 billion, Micro-C 2 billion, and Micro-C 3 billion data

Please refer to 1F for loosp calling. We called loops for different libraries (Hi-C 1 billion; Micro-C 1,2,3 billion data) at 1kb, 2kb, 5kb, 10kb, 25kb, 50kb resolutions respectively. The plot is showing the number of loops we got at different different resolutions among libraries.

### 2C. Fractions of loops that have different lengths found from Hi-C 1 billion, Micro-C 1 billion, 2 billion, and 3 billion data

We used a bash script ```category_loop_distance.sh``` to classify a specific loop data into different loop distance categoires. (<200kb, 200kb-400kb, 400kb-600kb, 600kb-800kb, 800kb-1Mb, > 1Mb) and summarize the fraction of each category. To use this script:

```
bash category_loop_distance.sh \
/PATH/TO/LOOPDATA \
OUTPUT_PATH
```

The default output folder path is your current working directory if you do not specify the output folder path. 
Outputs include a bunch of ```.bedpe``` files for loops falling into specific category, along with a ```summary.txt``` where the number and fraction of loops detected in different categoires are shown.

### 2D. Numbers of loops shared or unique among Hi-C 1 billion, Micro-C 1 billion, 2 billion, and 3 billion data

We used an Rscript ```Loop_comparison.R``` to run loop intersection analysis between two libraries, using a list of loop data as the input. 
See below code for usage instruction:

Load the function. To run it, make sure you have ```fuzzyjoin``` and ```dplyr``` libraries in your R environment:
```
source('Loop_comparison.R')
```

Then load the loop data. You should have at least two loop data to run this analysis. But let's use four as example here:
```
loop_1 = read.table(file = 'loop1.tsv')
loop_2 = read.table(file = 'loop2.tsv')
loop_3 = read.table(file = 'loop3.tsv')
loop_4 = read.table(file = 'loop4.tsv')
```

Create a NAMED list in your prefered order. Here, let's say the priority to assign overlap is loop_4 > loop_3 > loop_2 > loop_1, then:
```
loop_list <- list(
LOOP4 = loop_4,
LOOP3 = loop_3,
LOOP2 = loop_2,
LOOP1 = loop_1
)
```

Run the function using the list just created. The overlap will be calculated based on the priority same as the order in the list by default. The resolution argument is the resolution of your interaction data in bp.
```
loop_comparison(
loop_list = loop_list,
resolution = 5000,
output_path= './output'
)
```

Or if you want a different priority. For example, loop_1 > loop_2 > loop_3 > loop_4, then:
```
loop_comparison(
loop_list = loop_list,
priority = c(4,3,2,1),
resolution = 5000,
output_path= './output2'
)
```

The output files in the output folder include:
- All the original loop data added with comparison 0/1 matrix.
- For each loop data, output a subforder containning the shared/unique loop data
- Summary table for the loop comparison.

## Figure 3 Chromatin loops near structural variants

### 3A. Numbers of inter and intra-chromosomal structural variants identified from Hi-C and Micro-C data

We used hic_breakfinder (https://github.com/dixonlab/hic_breakfinder) to identify the structural variants. 

We used NeoLoopFinder (https://github.com/XiaoTaoWang/NeoLoopFinder, Wang et al, Nat Methods, 2021 PMID:34092790) to identify loops near the structural variants using Hi-C and Micro-C data. 

Please refere their github page for detailed instructions.

### 3B. Numbers of each category of structural variants identified from Hi-C and Micro-C data

Please refer to 3A. 

### 3C. Numbers of loops identified around the structural variants from Hi-C and Micro-C data

Please refer to 3A. 

### 3D. Numbers of neoloops that are shared or unique among Hi-C 1 billion, Micro-C 1 billion, 2 billion and 3 billion

Here is the loop-comparison analysis, we used the same comparison method as desribed in Figure 2D.

### 3E. An example heatmap of Micro-C data near the ARID1A gene that includes an inversion structural variant

We used a python script ```plotting_NEOheatmap.py``` to generate the heatmap together with multipe tracks of data. The necessary input to run this script include ```.cool``` with additional ```assemblies.txt``` and ```neo-loops.txt```, which are generate from NeoLoopFinder (Wang et al., Nat Methods, 2021, PMID:34092790). 

Below is an example code:
```
python plotting_NEOheatmap.py \
EXAMPLE.cool \
FromNeoLoopFinder.assemblies.txt \
FromNeoLoopFinder.neo-loops.txt \
./output/ \
Outputprefix
```
To plot with RNA-seq track, ```.bigwig``` file can be added in line 23 and 24 of the ```plotting_NEOheatmap.py```.

## Figure 4 Regulatory elements and nucleosome-depleted regions (NDRs) that are involved in loops

### 4A. Genome browser screenshots of ChIP-seq, NOMe-seq, Hi-C, Micro-C, and RefSeq Genes

This panel is generated using ```.bigwig``` files with [The Integrative Genomics Viewer (IGV)](https://software.broadinstitute.org/software/igv).

### 4B. Fractions of regulatory elements that intersect with loop anchors identified from Hi-C 1 billion, Micro-C 1 billion, 2 billion, and 3 billion data 

We used an Rscript ```check_REGelement.R``` to find out the regulatory elements involved in loops. To run this script, the loop data under specific resolution as well as the regulatory element data (in ```.bed``` format) is necessary. Check below code for example usage to identify enhancers involved in loop under 5kb resolution.
```
source('check_REGelement.R')
myloopdata <- read.table(file = '/PATH/TO/LoopData')
#E.g you want to find the enhancers involved in loops.
myEnhdata <- read.table(file = '/PATH/TO/EnhancerData')

check_REGelement(
loop_data = myloopdata, 
reg_data = myEnhdata,
reg_type = 'Enhancer', #To clarify, this option is just a character to help note the regulatory element you used.
resolution = 5000,
output_path='./output', 
prefix = 'outputprefix')
```
In the output, you will see two files ending with ```inLoop.bed``` or ```nonLoop.bed```. You can check ```summary.txt``` for the necessary stats.

### 4C. Numbers of loops belong to different loop categories defined by intersecting the loop anchors with different types of regulatory elements

We used an Rscript ```Loop_Rep_overlapAnalysis.R``` to perform the genomic regulatory element distribution analysis for loop data. See below instruction if you are interested to use this script.

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

### 4D. Comparison of number of promoter-enhancer loops identified from Hi-C 1 billion, Micro-C 1 billion, 2 billion and 3 billion data

Please refer to 4C for the analysis to identify promoter-enhancer loops.

### 4E. Significance of chromatin interaction (q-value identified by Mustache) for top 5 loop categories

Here we show an example Rscript ```plot_average_q_value.R``` which we used to generate the q-value violin plot of different chromatin interactions. Please note that the input ```.tsv``` files are generated from the ```Loop_Rep_overlapAnalysis.R``` script as described above.

## Figure 5 Promoter capture Micro-C data analysis 

### 5A. An overview of promoter capture Micro-C experimental procedure

### 5B. Chromatin interaction heatmaps of Micro-C and promoter capture Micro-C data near chr1q41 region

We used cooltools (https://github.com/open2c/cooltools) to generate the chromaitn interaction heatmap. We used the same script ```cooltools_heatmap.py``` as described in figure1.

### 5C. Significance of chromatin interaction for loops found in both promoter capture Micro-C and Micro-C (shared) and only one data 

We used the same script ```plot_average_q_value.R``` as described in Figure 4E to generate a violin plot.

### 5D. Fractions of active promoters that intersect with the loop anchors from Micro-C 1 billion, 2 billion, 3 billion data or promoter capture Micro-C data

Please refer to Figure 4B for the analysis.

### 5E. Numbers of loops and loop categories identified from promoter capture Micro-C data

We used the same script ```Loop_Rep_overlapAnalysis.R``` as described in Figure 4C.

## Figure 6 Nucleosome phasing and DNA methylation levels around regulatory elements involved in loops

### 6A. Average Micro-C signals around active promoters, enhancers, insulators, and NDRs without features that are in loop vs those that are not in loop

The signal plots are generated using DeepTools (Ramírez et al, Nucleic Acids Res, 2014 PMID:24799436). The input files used are ```.bigwig``` as well as the ```.bed``` files for regulatory elements (e.g promoter, enhancer, insulator, NDR). Please refer to GitHub page of DeepTools for detailed instructions: https://deeptools.readthedocs.io/en/develop/.

We used an Rscript ```t_test_4deeptools.R``` to perform t-test at the center of the signal between in-loop and not in-loop data. The ```RPKM.gz``` files were generated from DeepTools.

```
Rscript t_test_4deeptools.R /PATH/TO/inloop.RPKM.gz /PATH/TO/notinloop.RPKM.gz
```

### 6B. Average chromatin accessibility levels (%) of active promoters, enhancers, insulators, and NDRs without features that are in loop vs those that are not in loop 

This panel plots were generated using Bis-tools (Lay et al, Genome Res, 2015 PMID: 25747664). Please refer github page of Bis-tools for detailed instruction: https://github.com/dnaase/Bis-tools.
To make these graphs, the necessary input files, which are the ```HCG.bw```(DNA methylation) and ```GCH.bw```(chromatin accessibility), are generated from NOMe-seq data using Bis-SNP (Liu et al, Genome Biol, 2012 PMID: 22784381).

### 6C. Average DNA methylation levels of active promoters, enhancers, insulators, and NDRs without features that are in loop vs those that are not in loop 

Please refer to 6B.
