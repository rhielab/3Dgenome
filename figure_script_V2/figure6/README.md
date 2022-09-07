## Figure6 Nucleosome phasing and DNA methylation levels around regulatory elements involved in loops.

### 6A. Average Micro-C signals around active promoters, enhancers, insulators, and NDRs without features that are in loop vs those that are not in loop.

The signal plots are generated using DeepTools (Ramírez et al, 2014). The necessary input file is ```.bigwig``` as well as the ```.bed``` files for different regulatory elements (e.g promoter, enhancer, insulator, NDR). Please refer to https://deeptools.readthedocs.io/en/develop/ for more information. Below is our example code ro run DeepTools:

```
computeMatrix reference-point --referencePoint center \
-b 1000 -a 1000 \
-S input.bigwig \
-R RegulatoryElement.bed \
-p 16 \
-o output.RPKM.gz
```
The output ```RPKM.gz``` is necessary for ```plotHeatmap``` function
```
plotHeatmap \
-m output.RPKM.gz \
-o output.pdf
```
and also used to perform t-test. Here we provide an Rscript ```t_test_4deeptools.R``` to perform t-test at the center of the signal between in-loop and not in-loop data. To use:
```
Rscript t_test_4deeptools.R /PATH/TO/inloop.RPKM.gz /PATH/TO/notinloop.RPKM.gz
```
It will print out the t-test result.

### 6B. Average chromatin accessibility levels (%) of active promoters, enhancers, insulators, and NDRs without features that are in loop (black) vs those that are not in loop (orange) are shown. 6C. Average DNA methylation levels of active promoters, enhancers, insulators, and NDRs without features that are in loop (black) vs those that are not in loop (orange) are shown.

Panel C and D are generated using Bistools (Lay et al, 2015). Here we provide our example code to run Bistools. Basically, we used bistools to visualize NOMe-seq signal around interested sites (e.g. promoter, enhancer, insulator etc.) in density plot, average plot and heatmap.  Before the actual code to run Bisplot, some necessary variables need to be specified at first. Please noted that the ```HCG.bw``` and ```GCH.bw``` are generated from Bis-SNP ([Liu et al. 2012](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2012-13-7-r61)).
```
OUTPUTDIR=/PATH/TO/OUTPUTDIR
Prefix='OutPutprefix'
NAME='shortNAME' #short name like 'CTCF' or 'enhancer' or 'H3K4me3'
HCGBW='/PATH/TO/HCG.bw' #path to the HCG bigwig file from NOMEseq data
GCHBW='/PATH/TO/GCH.bw' #path to the GCG bigwig file from NOMEseq data
BED='/PATH/TO/GCH.bw' #a bed file to specifiy the genomic regions where you want to visualize the signal
```

To generate heatmap:
```
perl alignWigToBed.pl \
--density_bar \
--enrich_max 4.0 \
--result_dir $OUTPUTDIR \
--prefixs $Prefix \
--locs $L1 \
--category_names $NAME \
--sample_names $NAME \
--experiment_names Methylation --experiment_names Accessibility \
--rep_num_experiments 1 --rep_num_experiments 1 \
Samples.txt
```

To generate average plot:
```
perl alignWigToBed.pl \
HCG.bw GCH.bw \
--locs $BED \
--average \
--prefixs $Prefix \
--plot_x_axis_scale 1000 \
--data_matrix_scale 1200 \
--bin_size 20 \
--bin_size_align 1 \
--plot_x_axis_scale 1000 \
--result_dir $Dir \
--smooth \
--colors black --colors green \
--lengends $HCGBW --lengends $GCHBW
```

To genreate heatmap:
```
perl alignWigToBed.DKO_paper_version.pl \
--prefixs $Prefix \
--heatmap_with_reps \
--data_matrix_scale 1200 \
--bin_size 20 --bin_size_align 1 \
--plot_x_axis_scale 1000 \
--result_dir $OUTPUTDIR \
--locs $BED \
--category_names $NAME \
--experiment_names Methylation  \
--experiment_names Accessibility \
--experiment_names Methylation2  \
--experiment_names Accessibility2  \
--heatmap_col blue2yellow --heatmap_col white2darkgreen --heatmap_col blue2yellow --heatmap_col white2darkgreen \
--heatmap_regionToCluster_low -500 --heatmap_regionToCluster_high 500 \
--heatmap_cluster 2 \
--multiSampleClustering 4 \
Twosamples.txt
```

If you are interested to use BisTools, please check its github page for more usage instruction (https://github.com/dnaase/Bis-tools).
