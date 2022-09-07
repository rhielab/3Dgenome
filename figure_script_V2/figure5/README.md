## Figure5 Promoter capture Micro-C data analysis. 

### 5A. An overview of promoter capture Micro-C experimental procedure, including the promoter probe design scheme. Probes (green bar) with biotins (orange circle) are designed surrounding TSSs, and Micro-C reads are pulled down using the probes for promoter capture Micro-C.

### 5B. Chromatin interaction heatmaps of Micro-C and promoter capture Micro-C data near chr1q41 region at 2kb (top), 5kb (middle), and 10kb (bottom) resolutions.

We used cooltools to generate the chromaitn interaction heatmap. Please refer to [figure_script_V2/figure1](https://github.com/rhielab/3Dgenome/tree/main/figure_script_V2/figure1) for the usage.

### 5C. Significance of chromatin interaction (Chicago score (-log p-value), Mustache (q-value)) for loops found in both promoter capture Micro-C and Micro-C (shared) and only one data is plotted. A mean value in shown in red. A median value is shown in blue.

We used the same script ```plot_average_q_value.R``` as described in [figure_script_V2/figure4](https://github.com/rhielab/3Dgenome/tree/main/figure_script_V2/figure4) for violin plot generation.

### 5D. Fractions of active promoters that intersect with the loop anchors from Micro-C 1 billion, 2 billion, 3 billion data or promoter capture Micro-C data are shown (left). A fraction of active promoters that intersect with loop anchors from any datasets is shown in grey (In loop) while the one not in loop is shown in orange (Not in loop) (right).

We used an Rscript ```check_REGelement.R``` to find out the regulatory elements involved in loops. To run this script, the loop data under specific resolution as well as the regulatory element data (in ```.bed``` format) is necessary. Check below code for example usage to identify enhancers involved in loop under 5kb resolution.
```
source('check_REGelement.R')
myloopdata <- read.table(file = '/PATH/TO/LoopData')
#E.g you want to find the enhancers involved in loops.
myEnhdata <- read.table(file = '/PATH/TO/EnhancerData')

checkCGI_promoter(
loop_data = myloopdata, 
reg_data = myEnhdata,
reg_type = 'Enhancer', #To clarify, this option is just a character to help note the regulatory element you used.
resolution = 5000,
output_path='./output', 
prefix = 'outputprefix')
```
In the output, you will see two files ending with ```inLoop.bed``` or ```nonLoop.bed```. You can check ```summary.txt``` for the necessary stats.

### 5E. Numbers of loops and loop categories identified from promoter capture Micro-C data.

The analysis is bacially the stats from the regulatory elements distribution analysis among promoter capture Micro-C libraries. Please refer to [figure_script_v2/figure4](https://github.com/rhielab/3Dgenome/tree/main/figure_script_V2/figure4) for such analysis.
