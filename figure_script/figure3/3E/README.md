## Figure3 Chromatin loops near structural variants.

### 3E. An example heatmap of Micro-C data near the ARID1A gene that includes inversion structural variant is shown on the top. Under the heatmap, RNA-seq and RefSeq gene tracks are shown. Example neoloops newly gained due to the structural variants are circled in blue.

To generate the heatmap together with multipe tracks of data, here we provide a ```plotting_NEOheatmap.py``` script. The necessary input to run this script include ```.cool``` with additional ```assemblies.txt``` and ```neo-loops.txt``` which are generate from NeoLoopFinder (Wang et al., 2021). Please refert to [
/figure_script/figure3/3A.B.C](https://github.com/rhielab/3Dgenome/tree/main/figure_script/figure3/3A.B.C) for more details on how to generate these two files.

Below is an example code:
```
python plotting_NEOheatmap.py \
EXAMPLE.cool \
FromNeoLoopFinder.assemblies.txt \
FromNeoLoopFinder.neo-loops.txt \
./output/ \
Outputprefix
```

PLEASE NOTED that if you want to add tracks based on your need, the additional ```.bigwig``` files are needed and in this case you need to modify ```plotting_NEOheatmap.py``` by yourself (see line 23 and 24). Please refer to https://github.com/XiaoTaoWang/NeoLoopFinder for more details.
