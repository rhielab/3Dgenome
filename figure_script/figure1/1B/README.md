## Figure1 Comparison of Hi-C and Micro-C data.

### 1B. Chromatin interaction heatmaps of Hi-C and Micro-C data near chr7p14 region.

We used cooltools (https://github.com/open2c/cooltools) to generate the interaction heatmap. Here we provide a python script ```cooltools_heatmap.py``` to run cooltools. The required input to generate this plot is ```.mcool``` file. Below is the example code:

To run the script, please make sure that you have installed cooltools correctly. 
- 1st - mcool file
- 2nd - output path
- 3rd - prefix for the output
```
python cooltools_heatmap.py \
input/uni1945_2.5kb_1_billion.filt.mcool \
output/ \
20kb_chr7_31M-43M_Micro_C_1_Bil
```
