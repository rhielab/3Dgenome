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
and also used to perform t-test. We used an Rscript ```t_test_4deeptools.R``` to perform t-test at the center of the signal between in-loop and not in-loop data. To use:
```
Rscript t_test_4deeptools.R /PATH/TO/inloop.RPKM.gz /PATH/TO/notinloop.RPKM.gz
```
It will print out the t-test result.

### 6B. Average chromatin accessibility levels (%) of active promoters, enhancers, insulators, and NDRs without features that are in loop (black) vs those that are not in loop (orange) are shown. 6C. Average DNA methylation levels of active promoters, enhancers, insulators, and NDRs without features that are in loop (black) vs those that are not in loop (orange) are shown.

Panel C and D are generated using Bistools (Lay et al, 2015). Here we provide our example code to run Bistools. Basically, we used bistools to visualize NOMe-seq signal around interested sites (e.g. promoter, enhancer, insulator etc.) in density plot, average plot and heatmap.  To make these graphs, the necessary input files are the ```HCG.bw``` and ```GCH.bw```, which are generated from NOMe-seq data using Bis-SNP ([Liu et al. 2012](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2012-13-7-r61)), as well as a ```.bed``` file (e.g. Histong Modification or CTCF narrowPeak file) which is used to specify the regions that you want to plot the signals on.

To generate average plot:
```
perl alignWigToBed.pl \
HCG.bw GCH.bw \
--locs /PATH/TO/BED \
--average \
--prefixs output_prefix \
--plot_x_axis_scale 1000 \
--data_matrix_scale 1200 \
--bin_size 20 \
--bin_size_align 1 \
--plot_x_axis_scale 1000 \
--result_dir /path/to/outputfolder \
--smooth \
--colors black --colors green \
--lengends HCG.bw --lengends GCH.bw
```

To genreate heatmap, please noted that the ```shortName``` should be any features that describe your selected bed files (e.g Regulatory elements or Histone Marks):
```
perl alignWigToBed.DKO_paper_version.pl \
--prefixs output_prefix \
--heatmap_with_reps \
--data_matrix_scale 1200 \
--bin_size 20 --bin_size_align 1 \
--plot_x_axis_scale 1000 \
--result_dir /path/to/outputfolder \
--locs /PATH/TO/BED \
--category_names shortName \
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
