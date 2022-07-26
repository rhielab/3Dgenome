## **Figure-1 Comparison of Hi-C and Micro-C data.** 

**1B. Chromatin interaction heatmaps of Hi-C and Micro-C data near chr7p14 region.

We used cooltools (https://github.com/open2c/cooltools) to generate the interaction heatmap. Below is the example code to run cooltools:

```
#To run the script, require the cooltools conda environment. 
#1st - mcool file
#2nd - output path
#3rd - prefix for the output
python cooltools_heatmap.py \
input/uni1945_2.5kb_1_billion.filt.mcool \
output/ \
20kb_chr7_31M-43M_Micro_C_1_Bil
```
